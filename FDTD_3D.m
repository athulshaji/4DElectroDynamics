%% ----------------------------------------------------------------------
%   Program: FDTD_3D.m
%   Author:  Athul Shaji
%   Method: Finite Difference Time Domain - FDTD
%   Purpose: Wave scattering from 3D object(s) in the computation domain
%% ----------------------------------------------------------------------
clear all; clc; close all;
% type_sim=input('Simulation Type: 1)Monostatic 2) Bistatic? ');
type_sim=2;
% inc_scat_angles=input('incident/scattered angles in degrees [theta_i, phi_i, theta_s, phi_s]: ');
inc_scat_angles=[90,0,90,0];
%% Simulation Constants
c=2.99792458e8;   % Speed of Light
mu=4.0*pi*1.0e-7; % Free space permeability
eps=1.0/c^2/mu;   % Free space permittivity
epsr=1.0;           % Relative permittiviy
mur=1.0;            % Relative permeability2
eta=sqrt(mu*mur/eps/epsr); % Wave impedance
%% Simulation Enviroment Properties
pml_thickness=10; % PML thickness
pml_order=3; % PML order used for polynomial grading;
pml_rdesired=exp(-16); % PML desired maximum reflection error
pml_sigmafactor=1.3; % PML sigma factor used for grading sigma
pml_kappamax=7;  % PML maximum kappa used for grading kappa
pml_alphamin=0; %PML minimum alpha used for grading alpha
pml_alphamax=0.05; %PML maximum alpha used for grading alpha


if type_sim==1
    f_min=input('Enter Starting Frequency (f_begin)[GHz]:');
    f_max=input('Enter Ending Frequency (f_end)[GHz]:');
    f_step=input('Enter Frequency Step (f_delta)[GHz]:');
    freq=f_min:f_step:f_max;
elseif type_sim==2
    %     freq=input('Frequency(GHz)? ');
    freq=6;
    %     angle_step=input('Angle Step(Deg)? ');
    angle_step=1;
end
freq=freq*1e9;
incidentwave=[];  % Incident wave structure
incidentwave_ETheta = 1; % Incident wave theta component magnitude
incidentwave_EPhi = 0; % Incident wave phi component magnitude
incidentwave_inctheta = inc_scat_angles(1); % Incident wave theta angle
incidentwave_incphi = inc_scat_angles(2); % Incident wave phi angle

fmax=max(freq);  % Maximum observation frquency
lambdamin=c/fmax; % Minimum wavelength corresponds to maxium observation frequency

% N=input('Discretization resolution: (lambda/N); N=? '); % Number of cells per wavelength min
N=30;
delx=lambdamin/N; % Unit cell size in x direction
dely=delx; % Unit cell size in y direction
delz=delx; % Unit cell size in z direction

T=1000; % Total number of time steps

%% Simulation Materials
materials=[];   % Material Structure
% Free Space
materials(1).epsr=1; % Relative electrical permittivity
materials(1).mur=1; % Relative magnetic permittivity
materials(1).sigma=0; % Material electric conductivity
materials(1).sigmam=0; % Material magnetic conductivity

% PEC
materials(2).epsr=1; % Relative electrical permittivity
materials(2).mur=1; % Relative magnetic permittivity
materials(2).sigma=1e20; % Material electric conductivity
materials(2).sigmam=0; % Material magnetic conductivity

%% FDTD Domain Size

lx = -5e-2; % Rectangular Prism lower x coordinate
ly = -5e-2; % Rectangular Prism lower y coordinate
lz = -5e-2; % Rectangular Prism lower z coordinate

ux = 5e-2;  % Rectangular Prism upper x coordinate
uy = 5e-2;  % Rectangular Prism lower y coordinate
uz = 5e-2;  % Rectangular Prism lower z coordinate

freespace=5; % Number of free space cells added to the simulation domain

domain_lx = lx - delx * freespace -delx*pml_thickness; % Simulation Domain lower x coordinate with freespace and PML
domain_ly = ly - dely * freespace -dely*pml_thickness; % Simulation Domain lower y coordinate with freespace and PML
domain_lz = lz- delz * freespace - delz*pml_thickness; % Simulation Domain lower z coordinate with freespace and PML

domain_ux = ux+ delx * freespace + delx*pml_thickness; % Simulation Domain upper x coordinate with freespace and PML
domain_uy = uy+ dely * freespace + dely*pml_thickness; % Simulation Domain upper y coordinate with freespace and PML
domain_uz = uz+ delz * freespace + delz*pml_thickness; % Simulation Domain upper z coordinate with freespace and PML

[X,Y,Z]=meshgrid(domain_lx:delx:domain_ux,domain_ly:dely:domain_uy,domain_lz:delz:domain_uz);

nxtot = round((domain_ux - domain_lx)/delx);  % Total number of cells in x direction
nytot = round((domain_uy - domain_ly)/dely);  % Total number of cells in y direction
nztot = round((domain_uz - domain_lz)/delz);  % Total number of cells in z direction

domain_size_x = nxtot * delx; % Simulation size in x direction
domain_size_y = nytot * dely; % Simulation size in y direction
domain_size_z = nztot * delz; % Simulation size in z direction

domain_ux = domain_lx + domain_size_x; % Simulation Domain upper x coordinate with freespace and PML re-calculated because of rounding
domain_uy = domain_ly + domain_size_y; % Simulation Domain upper y coordinate with freespace and PML re-calculated because of rounding
domain_uz = domain_lz + domain_size_z; % Simulation Domain upper z coordinate with freespace and PML re-calculated because of rounding
%% Create Simulation Objects
simulationmaterials = ones(nxtot, nytot, nztot); % Simulation Material Array , initialized with free space
domain_centerx= zeros(nxtot,nytot,nztot); % Array for storing x direction cell center coordinates
domain_centery = zeros(nxtot,nytot,nztot); % Array for storing y direction cell center coordinates
domain_centerz = zeros(nxtot,nytot,nztot); % Array for storing z direction cell center coordinates

for ii = 1:nxtot
    domain_centerx(ii,:,:) = ...
        (ii - 0.5) * delx + domain_lx; % x direction cell center coordinates which is used for object creation
end
for jj = 1:nytot
    domain_centery(:,jj,:) = ...
        (jj - 0.5) * dely + domain_ly; % y direction cell center coordinates which is used for object creation
end
for kk = 1:nztot
    domain_centerz(:,:,kk) = ...
        (kk - 0.5) * delz + domain_lz; % z direction cell center coordinates which is used for object creation
end


% Assign Rectangular Prisms Materials

rectprismlx = round((lx - domain_lx)/delx) + 1;
rectprismly = round((ly - domain_ly)/dely) + 1;
rectprismlz = round((lz - domain_lz)/delz) + 1;

rectprismux = round((ux - domain_lx)/delx)+1;
rectprismuy = round((uy - domain_ly)/dely)+1;
rectprismuz = round((uz - domain_lz)/delz)+1;

simulationmaterials(rectprismlx:rectprismux-1, rectprismly:rectprismuy-1, rectprismlz:rectprismuz-1)=2;


%% Near Field Far Field - Frequency Domain Transformation Struct & Constitutive Parameters
ffield_distfromboundary=8; % Huygen Surface's Distance from the PML Absorbing Boundary

if type_sim==1
    ffield_nangles=1;    % Far field number of angles
elseif type_sim==2
    ffield_nangles=round(360/angle_step);    % Far field number of angles
end
ffield_plane='xy';  % Obsevation plane
ffield_theta=zeros(ffield_nangles,1); % Observation theta
ffield_phi=zeros(ffield_nangles,1); % Observation phi

ffield_dirTheta = zeros(size(freq,2),ffield_nangles); % Theta directivity
ffield_dir = zeros(size(freq,2),ffield_nangles);
ffield_dirPhi = zeros(size(freq,2),ffield_nangles); % Phi directivity


ffield_is=ffield_distfromboundary+1;  % Huygen Surface's x dimension beginning node                +--------+ (ie,je,ke)
ffield_js=ffield_distfromboundary+1;  % Huygen Surface's y dimension beginning node               /        /|
ffield_ks=ffield_distfromboundary+1;  % Huygen Surface's z dimension beginning node              +--------+ |
%                                                       |        | |
ffield_ie=nxtot-ffield_distfromboundary+1; % Huygen Surface's x dimesion termination node        |        | +
ffield_je=nytot-ffield_distfromboundary+1; % Huygen Surface's y dimesion termination node        |        | /
ffield_ke=nztot-ffield_distfromboundary+1; % Huygen Surface's z dimesion termination node        +--------+
%                                                (is,js,ks)

ffield_cx = 0.5*(ffield_is+ffield_ie); % Huygen Surface's center coordinate in x direction
ffield_cy = 0.5*(ffield_js+ffield_je); % Huygen Surface's center coordinate in y direction
ffield_cz = 0.5*(ffield_ks+ffield_ke); % Huygen Surface's center coordinate in z direction
%% Excitation Source Properties

source_x=round((nxtot+1)/2); % Excitation Source x coordinate
source_y=round((nytot+1)/2); % Excitation Source y coordinate
source_z=round(nztot/2); % Excitation Source z coordinate

tau=sqrt(2.3)/pi/fmax; % Parameter specifies the Gaussian waveform width in both time and frequency domain
t0=3*tau; % The amount of shift requires to avoid non-zero amplitudes at zero time.
J0=-1.0*eps;
%% Simulation Resolution Parameters

delt=delx/c/sqrt(3); % Increment in time (Cauran Stability Limit)

%% FDTD Update Coefficents Arrays
Ex=zeros(nxtot,nytot+1,nztot+1);  % Electric field x component array
Ey=zeros(nxtot+1,nytot,nztot+1);  % Electric field y component array
Ez=zeros(nxtot+1,nytot+1,nztot);  % Electric field z component array

Hx=zeros(nxtot+1,nytot,nztot);    % Magnetic field x component array
Hy=zeros(nxtot,nytot+1,nztot);    % Magnetic field y component array
Hz=zeros(nxtot,nytot,nztot+1);    % Magnetic field z component array

eps_x=ones(nxtot,nytot+1,nztot+1);   % Relative electrical permittivity array for x directiom
eps_y=ones(nxtot+1,nytot,nztot+1);   % Relative electrical permittivity array for y directiom
eps_z=ones(nxtot+1,nytot+1,nztot);   % Relative electrical permittivity array for z directiom

mu_x=ones(nxtot+1,nytot,nztot);      % Relative magnetic permeability array for x direction
mu_y=ones(nxtot,nytot+1,nztot);      % Relative magnetic permeability array for y direction
mu_z=ones(nxtot,nytot,nztot+1);      % Relative magnetic permeability array for z direction


sigma_x=zeros(nxtot,nytot+1,nztot+1); % Electrical conductivity array for x direction
sigma_y=zeros(nxtot+1,nytot,nztot+1); % Electrical conductivity array for y direction
sigma_z=zeros(nxtot+1,nytot+1,nztot); % Electrical conductivity array for z direction

sigmam_x=zeros(nxtot+1,nytot,nztot); % Magnetic conductivity array for x direction
sigmam_y=zeros(nxtot,nytot+1,nztot); % Magnetic conductivity array for y direction
sigmam_z=zeros(nxtot,nytot,nztot+1); % Magnetic conductivity array for z direction

Exincident_current = zeros(nxtot,nytot+1,nztot+1); % Array for storing Incident wave Ex component's current value
Eyincident_current = zeros(nxtot+1,nytot,nztot+1); % Array for storing Incident wave Ey component's current value
Ezincident_current = zeros(nxtot+1,nytot+1,nztot); % Array for storing Incident wave Ez component's current value

Exincident_past = zeros(nxtot,nytot+1,nztot+1); % Array for storing Incident wave Ex component's past value
Eyincident_past = zeros(nxtot+1,nytot,nztot+1); % Array for storing Incident wave Ey component's past value
Ezincident_past = zeros(nxtot+1,nytot+1,nztot); % Array for storing Incident wave Ez component's past value

Hxincident_current = zeros(nxtot+1,nytot,nztot); % Array for storing Incident wave Hx component's current value
Hyincident_current = zeros(nxtot,nytot+1,nztot); % Array for storing Incident wave Hy component's current value
Hzincident_current = zeros(nxtot,nytot,nztot+1); % Array for storing Incident wave Hz component's current value

Hxincident_past = zeros(nxtot+1,nytot,nztot); % Array for storing Incident wave Hx component's past value
Hyincident_past = zeros(nxtot,nytot+1,nztot); % Array for storing Incident wave Hy component's past value
Hzincident_past = zeros(nxtot,nytot,nztot+1); % Array for storing Incident wave Hz component's past value
%% Initialize Simulation Space Materials
% Four neighbor averaging used for calculating material properties
for ii = 1:size(materials,2)
    tempepsr(ii)   = materials(ii).epsr;
    tempmur(ii)    = materials(ii).mur;
    tempsigma(ii) = materials(ii).sigma;
    tempsigmam(ii) = materials(ii).sigmam;
end
tempmur(find(tempmur==0)) = 1e-20; % Avoid zero division
tempsigmam(find(tempsigmam==0)) = 1e-20; % Avoid zero division
%---------------------------------------------------------------------x components-----------------------------------------------------------------------
eps_x(1:nxtot,2:nytot,2:nztot)=0.25*(tempepsr(simulationmaterials(1:nxtot,2:nytot,2:nztot))+tempepsr(simulationmaterials(1:nxtot,1:nytot-1,2:nztot)) ...
    + tempepsr(simulationmaterials(1:nxtot,2:nytot,1:nztot-1))+ tempepsr(simulationmaterials(1:nxtot,1:nytot-1,1:nztot-1)));

mu_x(2:nxtot,1:nytot,1:nztot)=2*(tempmur(simulationmaterials(2:nxtot,1:nytot,1:nztot)).*tempmur(simulationmaterials(1:nxtot-1,1:nytot,1:nztot))) ...
    ./(tempmur(simulationmaterials(2:nxtot,1:nytot,1:nztot))+tempmur(simulationmaterials(1:nxtot-1,1:nytot,1:nztot)));

sigma_x(1:nxtot,2:nytot,2:nztot)=0.25 * (tempsigma(simulationmaterials(1:nxtot,2:nytot,2:nztot))+tempsigma(simulationmaterials(1:nxtot,1:nytot-1,2:nztot)) ...
    + tempsigma(simulationmaterials(1:nxtot,2:nytot,1:nztot-1))+tempsigma(simulationmaterials(1:nxtot,1:nytot-1,1:nztot-1)));

sigmam_x(2:nxtot,1:nytot,1:nztot)=2*(tempsigmam(simulationmaterials(2:nxtot,1:nytot,1:nztot)).*tempsigmam(simulationmaterials(1:nxtot-1,1:nytot,1:nztot))) ...
    ./(tempsigmam(simulationmaterials(2:nxtot,1:nytot,1:nztot))+tempsigmam(simulationmaterials(1:nxtot-1,1:nytot,1:nztot)));
%---------------------------------------------------------------------y components-----------------------------------------------------------------------
eps_y(2:nxtot,1:nytot,2:nztot)=0.25*(tempepsr(simulationmaterials(2:nxtot,1:nytot,2:nztot))+tempepsr(simulationmaterials(1:nxtot-1,1:nytot,2:nztot)) ...
    + tempepsr(simulationmaterials(2:nxtot,1:nytot,1:nztot-1))+tempepsr(simulationmaterials(1:nxtot-1,1:nytot,1:nztot-1)));

mu_y(1:nxtot,2:nytot,1:nztot)=2*(tempmur(simulationmaterials(1:nxtot,2:nytot,1:nztot)).*tempmur(simulationmaterials(1:nxtot,1:nytot-1,1:nztot))) ...
    ./(tempmur(simulationmaterials(1:nxtot,2:nytot,1:nztot))+tempmur(simulationmaterials(1:nxtot,1:nytot-1,1:nztot)));

sigma_y(2:nxtot,1:nytot,2:nztot)=0.25*(tempsigma(simulationmaterials(2:nxtot,1:nytot,2:nztot))+tempsigma(simulationmaterials(1:nxtot-1,1:nytot,2:nztot)) ...
    + tempsigma(simulationmaterials(2:nxtot,1:nytot,1:nztot-1))+tempsigma(simulationmaterials(1:nxtot-1,1:nytot,1:nztot-1)));

sigmam_y(1:nxtot,2:nytot,1:nztot)=2*(tempsigmam(simulationmaterials(1:nxtot,2:nytot,1:nztot)).*tempsigmam(simulationmaterials(1:nxtot,1:nytot-1,1:nztot))) ...
    ./(tempsigmam(simulationmaterials(1:nxtot,2:nytot,1:nztot))+tempsigmam(simulationmaterials(1:nxtot,1:nytot-1,1:nztot)));
%---------------------------------------------------------------------z components-----------------------------------------------------------------------
eps_z(2:nxtot,2:nytot,1:nztot)=0.25*(tempepsr(simulationmaterials(2:nxtot,2:nytot,1:nztot))+tempepsr(simulationmaterials(1:nxtot-1,2:nytot,1:nztot)) ...
    + tempepsr(simulationmaterials(2:nxtot,1:nytot-1,1:nztot))+tempepsr(simulationmaterials(1:nxtot-1,1:nytot-1,1:nztot)));

mu_z(1:nxtot,1:nytot,2:nztot)=2*(tempmur(simulationmaterials(1:nxtot,1:nytot,2:nztot)).*tempmur(simulationmaterials(1:nxtot,1:nytot,1:nztot-1))) ...
    ./(tempmur(simulationmaterials(1:nxtot,1:nytot,2:nztot))+tempmur(simulationmaterials(1:nxtot,1:nytot,1:nztot-1)));

sigma_z(2:nxtot,2:nytot,1:nztot)=0.25*(tempsigma(simulationmaterials(2:nxtot,2:nytot,1:nztot))+tempsigma(simulationmaterials(1:nxtot-1,2:nytot,1:nztot)) ...
    + tempsigma(simulationmaterials(2:nxtot,1:nytot-1,1:nztot))+tempsigma(simulationmaterials(1:nxtot-1,1:nytot-1,1:nztot)));

sigmam_z(1:nxtot,1:nytot,2:nztot)=2*(tempsigmam(simulationmaterials(1:nxtot,1:nytot,2:nztot)).*tempsigmam(simulationmaterials(1:nxtot,1:nytot,1:nztot-1))) ...
    ./(tempsigmam(simulationmaterials(1:nxtot,1:nytot,2:nztot))+tempsigmam(simulationmaterials(1:nxtot,1:nytot,1:nztot-1)));
%---------------------------------------------------------------------------------------------------------------------------------------------------------
%% Initialize Incident Field Update Coefficents
CEx_current= (2*(1-eps_x)*eps-delt*sigma_x)./(2*eps_x*eps+delt*sigma_x); % Coefficent used for Incident wave electric field x component current value
CEx_past=-(2*(1-eps_x)*eps+delt*sigma_x)./(2*eps_x*eps+delt*sigma_x); % Coefficent used for Incident wave electric field x component past value

CEy_current= (2*(1-eps_y)*eps-delt*sigma_y)./(2*eps_y*eps+delt*sigma_y); % Coefficent used for Incident wave electric field y component current value
CEy_past=-(2*(1-eps_y)*eps+delt*sigma_y)./(2*eps_y*eps+delt*sigma_y); % Coefficent used for Incident wave electric field y component past value

CEz_current= (2*(1-eps_z)*eps-delt*sigma_z)./(2*eps_z*eps+delt*sigma_z); % Coefficent used for Incident wave electric field z component current value
CEz_past=-(2*(1-eps_z)*eps+delt*sigma_z)./(2*eps_z*eps+delt*sigma_z); % Coefficent used for Incident wave electric field z component past value

CHx_current= (2*(1-mu_x)*mu-delt*sigmam_x)./(2*mu_x*mu+delt*sigmam_x); % Coefficent used for Incident wave magetic field x component current value
CHx_past=-(2*(1-mu_x)*mu+delt*sigmam_x)./(2*mu_x*mu+delt*sigmam_x); % Coefficent used for Incident wave magetic field x component past value

CHy_current= (2*(1-mu_y)*mu-delt*sigmam_y)./(2*mu_y*mu+delt*sigmam_y); % Coefficent used for Incident wave magetic field y component current value
CHy_past=-(2*(1-mu_y)*mu+delt*sigmam_y)./(2*mu_y*mu+delt*sigmam_y); % Coefficent used for Incident wave magetic field y component past value

CHz_current=(2*(1-mu_z)*mu-delt*sigmam_z)./(2*mu_z*mu+delt*sigmam_z); % Coefficent used for Incident wave magetic field z component current value
CHz_past=-(2*(1-mu_z)*mu+delt*sigmam_z)./(2*mu_z*mu+delt*sigmam_z); % Coefficent used for Incident wave magetic field z component past value

%% Initialize FDTD Updating Equations
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
C1Ex=(2*eps_x*eps-delt*sigma_x)./(2*eps_x*eps+delt*sigma_x); % First coefficent used to update Ex
C2Ex=(2*delt)/((2*eps_x*eps+delt*sigma_x)*dely); % Second coefficent used to update Ex
C3Ex=-(2*delt)/((2*eps_x*eps+delt*sigma_x)*delz); % Third coefficent used to update Ex

C1Ey=(2*eps_y*eps-delt*sigma_y)./(2*eps_y*eps+delt*sigma_y); % First coefficent used to update Ey
C2Ey=(2*delt)/((2*eps_y*eps+delt*sigma_y)*delz); % Second coefficent used to update Ey
C3Ey=-(2*delt)/((2*eps_y*eps+delt*sigma_y)*delx); % Third coefficent used to update Ey

C1Ez=(2*eps_z*eps-delt*sigma_z)./(2*eps_z*eps+delt*sigma_z); % First coefficent used to update Ez
C2Ez=(2*delt)/((2*eps_z*eps+delt*sigma_z)*delx); % Second coefficent used to update Ez
C3Ez=-(2*delt)/((2*eps_z*eps+delt*sigma_z)*dely); % Third coefficent used to update Ez

C1Hx=(2*mu_x*mu-delt*sigmam_x)./(2*mu_x*mu+delt*sigmam_x); % First coefficent used to update Hx
C2Hx=(2*delt)/((2*mu_x*mu+delt*sigmam_x)*delz); % Second coefficent used to update Hx
C3Hx=-(2*delt)/((2*mu_x*mu+delt*sigmam_x)*dely); % Third coefficent used to update Hx

C1Hy=(2*mu_y*mu-delt*sigmam_y)./(2*mu_y*mu+delt*sigmam_y); % First coefficent used to update Hy
C2Hy=(2*delt)/((2*mu_y*mu+delt*sigmam_y)*delx); % Second coefficent used to update Hy
C3Hy=-(2*delt)/((2*mu_y*mu+delt*sigmam_y)*delz); % Third coefficent used to update Hy

C1Hz=(2*mu_z*mu-delt*sigmam_z)./(2*mu_z*mu+delt*sigmam_z); % First coefficent used to update Hz
C2Hz=(2*delt)/((2*mu_z*mu+delt*sigmam_z)*dely); % Second coefficent used to update Ez
C3Hz=-(2*delt)/((2*mu_z*mu+delt*sigmam_z)*delx); % Third coefficent used to update Ez
%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%% Incident Field Parameters
inctheta = incidentwave_inctheta*pi/180;
incphi = incidentwave_incphi*pi/180;
ETheta = incidentwave_ETheta;
EPhi = incidentwave_EPhi;

Exincamp = ETheta * cos(inctheta) * cos(incphi)- EPhi * sin(incphi);
Eyincamp = ETheta * cos(inctheta) * sin(incphi)+ EPhi * cos(incphi);
Ezincamp = -ETheta * sin(inctheta);

Hxincamp = (-1/eta)*(EPhi * cos(inctheta)*cos(incphi) + ETheta * sin(incphi));
Hyincamp = (-1/eta)*(EPhi * cos(inctheta)*sin(incphi) - ETheta * cos(incphi));
Hzincamp = (1/eta)*(EPhi * sin(inctheta));

cordx = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing x coordinates
cordy = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing y coordinates
cordz = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing z coordinates
for ii= 1:nxtot+1
    cordx(ii,:,:) =  domain_lx + (ii - 1) * delx;
end
for jj = 1:nytot+1
    cordy(:,jj,:) = domain_ly + (jj - 1) * dely;
end
for kk = 1:nztot+1
    cordz(:,:,kk) = domain_lz+(kk - 1) * delz;
end


d0 =[domain_lx domain_ly domain_lz; domain_lx domain_ly domain_uz; domain_lx domain_uy domain_lz; domain_lx domain_uy domain_uz;
    domain_ux domain_ly domain_lz; domain_ux domain_ly domain_uz;  domain_ux domain_uy domain_lz; domain_ux domain_uy domain_uz;];

k_x =  sin(inctheta)*cos(incphi); % k vector x component
k_y =  sin(inctheta)*sin(incphi); % k vector y component
k_z =  cos(inctheta); % k vector z component

kr0 = k_x * d0(:,1)+k_y * d0(:,2)+k_z * d0(:,3);
spatialshift= min(kr0)/c; % Calculate Spatial Shift

%--------------------------------------------------------- calculate k.r for every field component--------------------------------------------------------
krEx = -spatialshift+((cordx(1:nxtot,1:nytot+1,1:nztot+1)+delx/2) * k_x + cordy(1:nxtot,1:nytot+1,1:nztot+1) * k_y + cordz(1:nxtot,1:nytot+1,1:nztot+1)...
    * k_z)/c;
krEy = -spatialshift+(cordx(1:nxtot+1,1:nytot,1:nztot+1) * k_x + (cordy(1:nxtot+1,1:nytot,1:nztot+1)+dely/2) * k_y + cordz(1:nxtot+1,1:nytot,1:nztot+1)...
    * k_z)/c;
krEz = -spatialshift+(cordx(1:nxtot+1,1:nytot+1,1:nztot) * k_x + cordy(1:nxtot+1,1:nytot+1,1:nztot) * k_y + (cordz(1:nxtot+1,1:nytot+1,1:nztot)+delz/2)...
    * k_z)/c;
krHx = -spatialshift+(cordx(1:nxtot+1,1:nytot,1:nztot) * k_x + (cordy(1:nxtot+1,1:nytot,1:nztot)+dely/2) * k_y + (cordz(1:nxtot+1,1:nytot,1:nztot)...
    +delz/2) * k_z)/c;
krHy = -spatialshift+((cordx(1:nxtot,1:nytot+1,1:nztot)+delx/2) * k_x + cordy(1:nxtot,1:nytot+1,1:nztot) * k_y + (cordz(1:nxtot,1:nytot+1,1:nztot)...
    +delz/2) * k_z)/c;
krHz = -spatialshift+((cordx(1:nxtot,1:nytot,1:nztot+1)+delx/2) * k_x + (cordy(1:nxtot,1:nytot,1:nztot+1)+dely/2) * k_y + cordz(1:nxtot,1:nytot,1:nztot+1)...
    * k_z)/c;

clear cordx cordy cordz;

%% Initialize CPML Parameters
sigmaopt=(pml_order+1)/(150*pi*delx); % Optimal sigma used for grading conductivity in x direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in x direction

% xleft face sigma(x) and sigmam(x), kappa(x) , kappam(x), alpha(x) and alpham(x) grading

% The sigma(x),kappa(x),alpha(x) are used to update field components Hy and Hz , sigmam(x),kappam(x) and alpham(x)
% are used to update field components Ez and Ey.

% The positions of Hy and Hz in x direction has an offset of (i+1/2) , on the other hand Ez and Ey has an offset of (i)


%----------------------------------------------------------------------xleft region-------------------------------------------------------------------------
Psi_ey_hz_xl=zeros(pml_thickness,nytot,nztot+1);
Psi_ez_hy_xl=zeros(pml_thickness,nytot+1,nztot);


Psi_hy_ez_xl=zeros(pml_thickness,nytot+1,nztot);
Psi_hz_ey_xl=zeros(pml_thickness,nytot,nztot+1);

xm=zeros(1,pml_thickness);  % Magnetic field components position
xe=zeros(1,pml_thickness);  % Electric field components position

sigmax=zeros(1,pml_thickness);  % Electrical conductivity profile in x direction
sigmamx=zeros(1,pml_thickness);  % Magnetic conductivity profile in x direction

kappax=zeros(1,pml_thickness);  % Parameter kappa_e profile in x direction
kappamx=zeros(1,pml_thickness); % Parameter kappa_m in x direction

alphax=zeros(1,pml_thickness);  % Parameter alpha_e profile in x direction
alphamx=zeros(1,pml_thickness);  % Parameter alpham_e profile in x direction

for ii=1:pml_thickness
    xm(ii)=(pml_thickness-ii+0.75)/pml_thickness;
    xe(ii)=(pml_thickness-ii+0.25)/pml_thickness;
    
    sigmax(ii)=sigmamax*xe(ii)^pml_order;
    sigmamx(ii)=(mu/eps)*sigmamax*xm(ii)^pml_order;
    
    
    kappax(ii)=1+(pml_kappamax-1)*xe(ii)^pml_order;
    kappamx(ii)=1+(pml_kappamax-1)*xm(ii)^pml_order;
    
    alphax(ii)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-xe(ii));
    alphamx(ii)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-xm(ii)));
end

bxl = exp((-delt/eps)*((sigmax./kappax)+ alphax));   % The parameter b which is used to calculate Psi auxilary field
bxml= exp((-delt/mu)*((sigmamx./kappamx)+ alphamx)); % The parameter b_m which is used to calculate Psi auxilary field
axl= (1/delx)*(bxl-1.0).* sigmax./(kappax.*(sigmax+kappax.*alphax)); % The parameter a which is used to calculate Psi auxilary field
axml= (1/delx)*(bxml-1.0).* sigmamx./(kappamx.*(sigmamx+kappamx.*alphamx)); % The parameter a which is used to calculate Psi auxilary field


CPsi_ey_hz_xl=C3Ey(2:pml_thickness+1,:,:)*delx;
CPsi_ez_hy_xl=C2Ez(2:pml_thickness+1,:,:)*delx;

CPsi_hy_ez_xl=C2Hy(1:pml_thickness,:,:)*delx;
CPsi_hz_ey_xl=C3Hz(1:pml_thickness,:,:)*delx;


for ii = 1: pml_thickness
    C3Ey(ii+1,:,:) = C3Ey(ii+1,:,:)/kappax(ii);
    C2Ez(ii+1,:,:) = C2Ez(ii+1,:,:)/kappax(ii);
    C2Hy(ii,:,:) = C2Hy(ii,:,:)/kappamx(ii);
    C3Hz(ii,:,:) = C3Hz(ii,:,:)/kappamx(ii);
end
%-----------------------------------------------------------------End of xleft region-----------------------------------------------------------------------

%----------------------------------------------------------------------xright region-------------------------------------------------------------------------
Psi_ey_hz_xr=zeros(pml_thickness,nytot,nztot+1);
Psi_ez_hy_xr=zeros(pml_thickness,nytot+1,nztot);

Psi_hy_ez_xr=zeros(pml_thickness,nytot+1,nztot);
Psi_hz_ey_xr=zeros(pml_thickness,nytot,nztot+1);


xm=zeros(1,pml_thickness);  % Magnetic field components position
xe=zeros(1,pml_thickness);  % Electric field components position

sigmax=zeros(1,pml_thickness);  % Electrical conductivity profile in x direction
sigmamx=zeros(1,pml_thickness);  % Magnetic conductivity profile in x direction

kappax=zeros(1,pml_thickness);  % Parameter kappa_e profile in x direction
kappamx=zeros(1,pml_thickness); % Parameter kappa_m in x direction

alphax=zeros(1,pml_thickness);  % Parameter alpha_e profile in x direction
alphamx=zeros(1,pml_thickness);  % Parameter alpham_e profile in x direction

for ii=1:pml_thickness
    xe(ii)=(ii-0.75)/pml_thickness;
    xm(ii)=(ii-0.25)/pml_thickness;
    
    sigmax(ii)=sigmamax*xe(ii)^pml_order;
    sigmamx(ii)=(mu/eps)*sigmamax*xm(ii)^pml_order;
    
    
    kappax(ii)=1+(pml_kappamax-1)*xe(ii)^pml_order;
    kappamx(ii)=1+(pml_kappamax-1)*xm(ii)^pml_order;
    
    alphax(ii)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-xe(ii));
    alphamx(ii)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-xm(ii)));
end

bxr = exp((-delt/eps)*((sigmax./kappax)+ alphax));   % The parameter b which is used to calculate Psi auxilary field
bxmr= exp((-delt/mu)*((sigmamx./kappamx)+ alphamx)); % The parameter b_m which is used to calculate Psi auxilary field
axr= (1/delx)*(bxr-1.0).* sigmax./(kappax.*(sigmax+kappax.*alphax)); % The parameter a which is used to calculate Psi auxilary field
axmr= (1/delx)*(bxmr-1.0).* sigmamx./(kappamx.*(sigmamx+kappamx.*alphamx)); % The parameter a which is used to calculate Psi auxilary field



CPsi_ey_hz_xr=C3Ey(nxtot+1-pml_thickness:nxtot,:,:)*delx;
CPsi_ez_hy_xr=C2Ez(nxtot+1-pml_thickness:nxtot,:,:)*delx;

CPsi_hy_ez_xr=C2Hy(nxtot+1-pml_thickness:nxtot,:,:)*delx;
CPsi_hz_ey_xr=C3Hz(nxtot+1-pml_thickness:nxtot,:,:)*delx;


for ii = 1: pml_thickness
    C3Ey(nxtot-pml_thickness+ii,:,:) = C3Ey(nxtot-pml_thickness+ii,:,:)/kappax(ii);
    C2Ez(nxtot-pml_thickness+ii,:,:) = C2Ez(nxtot-pml_thickness+ii,:,:)/kappax(ii);
    C2Hy(nxtot-pml_thickness+ii,:,:) = C2Hy(nxtot-pml_thickness+ii,:,:)/kappamx(ii);
    C3Hz(nxtot-pml_thickness+ii,:,:) = C3Hz(nxtot-pml_thickness+ii,:,:)/kappamx(ii);
end
%-----------------------------------------------------------------End of xright region-----------------------------------------------------------------------

%--------------------------------------------------------------------yleft region----------------------------------------------------------------------------

sigmaopt=(pml_order+1)/(150*pi*dely); % Optimal sigma used for grading conductivity in y direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in y direction

Psi_ex_hz_yl = zeros(nxtot,pml_thickness,nztot+1);
Psi_ez_hx_yl = zeros(nxtot+1,pml_thickness,nztot);

Psi_hx_ez_yl = zeros(nxtot+1,pml_thickness,nztot);
Psi_hz_ex_yl = zeros(nxtot,pml_thickness,nztot+1);

ym=zeros(1,pml_thickness);  % Magnetic field components position
ye=zeros(1,pml_thickness);  % Electric field components position

sigmay=zeros(1,pml_thickness);  % Electrical conductivity profile in y direction
sigmamy=zeros(1,pml_thickness);  % Magnetic conductivity profile in y direction

kappay=zeros(1,pml_thickness);  % Parameter kappa_e profile in y direction
kappamy=zeros(1,pml_thickness); % Parameter kappa_m in y direction

alphay=zeros(1,pml_thickness);  % Parameter alpha_e profile in y direction
alphamy=zeros(1,pml_thickness);  % Parameter alpham_e profile in y direction

for jj=1:pml_thickness
    ym(jj)=(pml_thickness-jj+0.75)/pml_thickness; %Magnetic field components position calculation
    ye(jj)=(pml_thickness-jj+0.25)/pml_thickness; %Electric field components position calculation
    
    sigmay(jj)=sigmamax*ye(jj)^pml_order;
    sigmamy(jj)=(mu/eps)*sigmamax*ym(jj)^pml_order;
    
    
    kappay(jj)=1+(pml_kappamax-1)*ye(jj)^pml_order;
    kappamy(jj)=1+(pml_kappamax-1)*ym(jj)^pml_order;
    
    alphay(jj)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ye(jj));
    alphamy(jj)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ym(jj)));
end

byl = exp((-delt/eps)*((sigmay./kappay)+ alphay));   % The parameter b which is used to calculate Psi auxilary field
byml= exp((-delt/mu)*((sigmamy./kappamy)+ alphamy)); % The parameter b_m which is used to calculate Psi auxilary field
ayl= (1/dely)*(byl-1.0).* sigmay./(kappay.*(sigmay+kappay.*alphay)); % The parameter a which is used to calculate Psi auxilary field
ayml= (1/dely)*(byml-1.0).* sigmamy./(kappamy.*(sigmamy+kappamy.*alphamy)); % The parameter a which is used to calculate Psi auxilary field


CPsi_ex_hz_yl = C2Ex(:,2:pml_thickness+1,:)*dely;
CPsi_ez_hx_yl = C3Ez(:,2:pml_thickness+1,:)*dely;

CPsi_hx_ez_yl = C3Hx(:,1:pml_thickness,:)*dely;
CPsi_hz_ex_yl = C2Hz(:,1:pml_thickness,:)*dely;



for jj = 1: pml_thickness
    C2Ex(:,jj+1,:) = C2Ex(:,jj+1,:)/kappay(jj);
    C3Ez(:,jj+1,:) = C3Ez(:,jj+1,:)/kappay(jj);
    C3Hx(:,jj,:) = C3Hx(:,jj,:)/kappamy(jj);
    C2Hz(:,jj,:) = C2Hz(:,jj,:)/kappamy(jj);
end
%-----------------------------------------------------------------End of yleft region-----------------------------------------------------------------------

%--------------------------------------------------------------------yright region-------------------------------------------------------------------------

sigmaopt=(pml_order+1)/(150*pi*dely); % Optimal sigma used for grading conductivity in y direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in y direction

Psi_ex_hz_yr = zeros(nxtot,pml_thickness,nztot+1);
Psi_ez_hx_yr = zeros(nxtot+1,pml_thickness,nztot);

Psi_hx_ez_yr = zeros(nxtot+1,pml_thickness,nztot);
Psi_hz_ex_yr = zeros(nxtot,pml_thickness,nztot+1);

ym=zeros(1,pml_thickness);  % Magnetic field components position
ye=zeros(1,pml_thickness);  % Electric field components position

sigmay=zeros(1,pml_thickness);  % Electrical conductivity profile in y direction
sigmamy=zeros(1,pml_thickness);  % Magnetic conductivity profile in y direction

kappay=zeros(1,pml_thickness);  % Parameter kappa_e profile in y direction
kappamy=zeros(1,pml_thickness); % Parameter kappa_m in y direction

alphay=zeros(1,pml_thickness);  % Parameter alpha_e profile in y direction
alphamy=zeros(1,pml_thickness);  % Parameter alpham_e profile in y direction

for jj=1:pml_thickness
    ye(jj)=(jj-0.75)/pml_thickness;
    ym(jj)=(jj-0.25)/pml_thickness;
    
    sigmay(jj)=sigmamax*ye(jj)^pml_order;
    sigmamy(jj)=(mu/eps)*sigmamax*ym(jj)^pml_order;
    
    
    kappay(jj)=1+(pml_kappamax-1)*ye(jj)^pml_order;
    kappamy(jj)=1+(pml_kappamax-1)*ym(jj)^pml_order;
    
    alphay(jj)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ye(jj));
    alphamy(jj)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ym(jj)));
end

byr = exp((-delt/eps)*((sigmay./kappay)+ alphay));   % The parameter b which is used to calculate Psi auxilary field
bymr= exp((-delt/mu)*((sigmamy./kappamy)+ alphamy)); % The parameter b_m which is used to calculate Psi auxilary field
ayr= (1/dely)*(byr-1.0).* sigmay./(kappay.*(sigmay+kappay.*alphay)); % The parameter a which is used to calculate Psi auxilary field
aymr= (1/dely)*(bymr-1.0).* sigmamy./(kappamy.*(sigmamy+kappamy.*alphamy)); % The parameter a which is used to calculate Psi auxilary field


CPsi_ex_hz_yr = C2Ex(:,nytot+1-pml_thickness:nytot,:)*dely;
CPsi_ez_hx_yr = C3Ez(:,nytot+1-pml_thickness:nytot,:)*dely;

CPsi_hx_ez_yr = C3Hx(:,nytot+1-pml_thickness:nytot,:)*dely;
CPsi_hz_ex_yr = C2Hz(:,nytot+1-pml_thickness:nytot,:)*dely;

for jj = 1: pml_thickness
    C2Ex(:,nytot-pml_thickness+jj,:) = C2Ex(:,nytot-pml_thickness+jj,:)/kappay(jj);
    C3Ez(:,nytot-pml_thickness+jj,:) = C3Ez(:,nytot-pml_thickness+jj,:)/kappay(jj);
    C3Hx(:,nytot-pml_thickness+jj,:) = C3Hx(:,nytot-pml_thickness+jj,:)/kappamy(jj);
    C2Hz(:,nytot-pml_thickness+jj,:) = C2Hz(:,nytot-pml_thickness+jj,:)/kappamy(jj);
end
%-----------------------------------------------------------------End of yright region-----------------------------------------------------------------------

%--------------------------------------------------------------------zleft region----------------------------------------------------------------------------
sigmaopt=(pml_order+1)/(150*pi*delz); % Optimal sigma used for grading conductivity in z direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in z direction

Psi_ex_hy_zl = zeros(nxtot,nytot+1,pml_thickness);
Psi_ey_hx_zl = zeros(nxtot+1,nytot,pml_thickness);

Psi_hx_ey_zl = zeros(nxtot+1,nytot,pml_thickness);
Psi_hy_ex_zl = zeros(nxtot,nytot+1,pml_thickness);

zm=zeros(1,pml_thickness);  % Magnetic field components position
ze=zeros(1,pml_thickness);  % Electric field components position

sigmaz=zeros(1,pml_thickness);  % Electrical conductivity profile in z direction
sigmamz=zeros(1,pml_thickness);  % Magnetic conductivity profile in z direction

kappaz=zeros(1,pml_thickness);  % Parameter kappa_e profile in z direction
kappamz=zeros(1,pml_thickness); % Parameter kappa_m in z direction

alphaz=zeros(1,pml_thickness);  % Parameter alpha_e profile in z direction
alphamz=zeros(1,pml_thickness);  % Parameter alpham_e profile in z direction

for kk=1:pml_thickness
    zm(kk)=(pml_thickness-kk+0.75)/pml_thickness;
    ze(kk)=(pml_thickness-kk+0.25)/pml_thickness;
    
    sigmaz(kk)=sigmamax*ze(kk)^pml_order;
    sigmamz(kk)=(mu/eps)*sigmamax*zm(kk)^pml_order;
    
    kappaz(kk)=1+(pml_kappamax-1)*ze(kk)^pml_order;
    kappamz(kk)=1+(pml_kappamax-1)*zm(kk)^pml_order;
    
    alphaz(kk)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ze(kk));
    alphamz(kk)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-zm(kk)));
end

bzl = exp((-delt/eps)*((sigmaz./kappaz)+ alphaz));   % The parameter b which is used to calculate Psi auxilary field
bzml= exp((-delt/mu)*((sigmamz./kappamz)+ alphamz)); % The parameter b_m which is used to calculate Psi auxilary field
azl= (1/delz)*(bzl-1.0).* sigmaz./(kappaz.*(sigmaz+kappaz.*alphaz)); % The parameter a which is used to calculate Psi auxilary field
azml= (1/delz)*(bzml-1.0).* sigmamz./(kappamz.*(sigmamz+kappamz.*alphamz)); % The parameter a which is used to calculate Psi auxilary field


CPsi_ex_hy_zl = C3Ex(:,:,2:pml_thickness+1)*delz;
CPsi_ey_hx_zl = C2Ey(:,:,2:pml_thickness+1)*delz;

CPsi_hx_ey_zl = C2Hx(:,:,1:pml_thickness)*delz;
CPsi_hy_ex_zl = C3Hy(:,:,1:pml_thickness)*delz;

for kk = 1: pml_thickness
    C3Ex(:,:,kk+1) = C3Ex(:,:,kk+1)/kappaz(kk);
    C2Ey(:,:,kk+1) = C2Ey(:,:,kk+1)/kappaz(kk);
    
    C2Hx(:,:,kk) = C2Hx(:,:,kk)/kappamz(kk);
    C3Hy(:,:,kk) = C3Hy(:,:,kk)/kappamz(kk);
end
%-----------------------------------------------------------------End of zleft region-----------------------------------------------------------------------

%--------------------------------------------------------------------zright region----------------------------------------------------------------------------
sigmaopt=(pml_order+1)/(150*pi*delz); % Optimal sigma used for grading conductivity in z direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in z direction

Psi_ex_hy_zr = zeros(nxtot,nytot+1,pml_thickness);
Psi_ey_hx_zr = zeros(nxtot+1,nytot,pml_thickness);

Psi_hx_ey_zr = zeros(nxtot+1,nytot,pml_thickness);
Psi_hy_ex_zr = zeros(nxtot,nytot+1,pml_thickness);

zm=zeros(1,pml_thickness);  % Magnetic field components position
ze=zeros(1,pml_thickness);  % Electric field components position

sigmaz=zeros(1,pml_thickness);  % Electrical conductivity profile in z direction
sigmamz=zeros(1,pml_thickness);  % Magnetic conductivity profile in z direction

kappaz=zeros(1,pml_thickness);  % Parameter kappa_e profile in z direction
kappamz=zeros(1,pml_thickness); % Parameter kappa_m in z direction

alphaz=zeros(1,pml_thickness);  % Parameter alpha_e profile in z direction
alphamz=zeros(1,pml_thickness);  % Parameter alpham_e profile in z direction

for kk=1:pml_thickness
    ze(kk)=(kk-0.75)/pml_thickness;
    zm(kk)=(kk-0.25)/pml_thickness;
    
    sigmaz(kk)=sigmamax*ze(kk)^pml_order;
    sigmamz(kk)=(mu/eps)*sigmamax*zm(kk)^pml_order;
    
    kappaz(kk)=1+(pml_kappamax-1)*ze(kk)^pml_order;
    kappamz(kk)=1+(pml_kappamax-1)*zm(kk)^pml_order;
    
    alphaz(kk)=pml_alphamin+(pml_alphamax-pml_alphamin)*(1-ze(kk));
    alphamz(kk)=(mu/eps)*(pml_alphamin+(pml_alphamax-pml_alphamin)*(1-zm(kk)));
end

bzr = exp((-delt/eps)*((sigmaz./kappaz)+ alphaz));   % The parameter b which is used to calculate Psi auxilary field
bzmr= exp((-delt/mu)*((sigmamz./kappamz)+ alphamz)); % The parameter b_m which is used to calculate Psi auxilary field
azr= (1/delz)*(bzr-1.0).* sigmaz./(kappaz.*(sigmaz+kappaz.*alphaz)); % The parameter a which is used to calculate Psi auxilary field
azmr= (1/delz)*(bzmr-1.0).* sigmamz./(kappamz.*(sigmamz+kappamz.*alphamz)); % The parameter a which is used to calculate Psi auxilary field


CPsi_ex_hy_zr = C3Ex(:,:,nztot+1-pml_thickness:nztot)*delz;
CPsi_ey_hx_zr = C2Ey(:,:,nztot+1-pml_thickness:nztot)*delz;

CPsi_hx_ey_zr = C2Hx(:,:,nztot+1-pml_thickness:nztot)*delz;
CPsi_hy_ex_zr = C3Hy(:,:,nztot+1-pml_thickness:nztot)*delz;

for kk = 1: pml_thickness
    C3Ex(:,:,nztot-pml_thickness+kk) = C3Ex(:,:,nztot-pml_thickness+kk)/kappaz(kk);
    C2Ey(:,:,nztot-pml_thickness+kk) = C2Ey(:,:,nztot-pml_thickness+kk)/kappaz(kk);
    
    C2Hx(:,:,nztot-pml_thickness+kk) = C2Hx(:,:,nztot-pml_thickness+kk)/kappamz(kk);
    C3Hy(:,:,nztot-pml_thickness+kk) = C3Hy(:,:,nztot-pml_thickness+kk)/kappamz(kk);
end
%-----------------------------------------------------------------End of zright region-----------------------------------------------------------------------

%---------------------------------------------------------------End of CPML Initialization -----------------------------------------------------------------

%% --------------------------------------------Initialize time domain electric andmagnetic current arrays-----------------------------------------
% For +x face we have Jy,Jz, My,Mz according to J = x x H ; and M = -x x E;
tdjyxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current y component array calculated at +x surface
tdjzxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at +x surface
tdmyxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current y component array calculated at +x surface
tdmzxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at +x surface

% For -x face we have Jy,Jz, My,Mz according to J = -x x H ; and M = -(-x) x E;
tdjyxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current y component array calculated at -x surface
tdjzxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at -x surface
tdmyxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current y component array calculated at -x surface
tdmzxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at -x surface

%------------------------------------------------------------------------------------------------------------------------------------------

% For +y face we have Jx,Jz,Mx,Mz according to J = y x H ; and M = -y x E;
tdjxyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current x component array calculated at +y surface
tdjzyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at +y surface
tdmxyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current x component array calculated at +y surface
tdmzyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at +y surface

% For -y face we have Jx,Jz,Mx,Mz according to J = -y x H ; and M = -(-y) x E;
tdjxyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current x component array calculated at -y surface
tdjzyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at -y surface
tdmxyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current x component array calculated at -y surface
tdmzyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at -y surface

%------------------------------------------------------------------------------------------------------------------------------------------

% For +z face we have Jx,Jy, Mx,My according to J = z x H ; and M = -z x E;
tdjxzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current x component array calculated at +z surface
tdjyzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current y component array calculated at +z surface
tdmxzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current x component array calculated at +z surface
tdmyzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current y component array calculated at +z surface

% For -z face we have Jx,Jy, Mx,My according to J = -z x H ; and M = -(-z) x E;
tdjxzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current x component array calculated at -z surface
tdjyzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current y component array calculated at -z surface
tdmxzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current x component array calculated at -z surface
tdmyzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current y component array calculated at -z surface

%------------------------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------End of time domain electric andmagnetic current arrays-----------------------------------------

%--------------------------------------------Frequency domain electric and magnetic current arrays-----------------------------------------

% The first dimension in these array reserved for DFT frequencies

%------------------------------------------------------------------------------------------------------------------------------------------

% For +x face we have Jy,Jz, My,Mz according to J = x x H ; and M = -x x E;
fdjyxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current y component array calculated at +x surface
fdjzxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at +x surface
fdmyxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current y component array calculated at +x surface
fdmzxp=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at +x surface

% For -x face we have Jy,Jz, My,Mz according to J = -x x H ; and M = -(-x) x E;
fdjyxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current y component array calculated at -x surface
fdjzxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at -x surface
fdmyxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current y component array calculated at -x surface
fdmzxn=zeros(size(freq,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at -x surface

%----------------------------------------------------------------------------------------------------------------------------------------

% For +y face we have Jx,Jz,Mx,Mz according to J = y x H ; and M = -y x E;
fdjxyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current x component array calculated at +y surface
fdjzyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at +y surface
fdmxyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current x component array calculated at +y surface
fdmzyp=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at +y surface

% For -y face we have Jx,Jz,Mx,Mz according to J = -y x H ; and M = -(-y) x E;
fdjxyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current x component array calculated at -y surface
fdjzyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at -y surface
fdmxyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current x component array calculated at -y surface
fdmzyn=zeros(size(freq,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at -y surface

%--------------------------------------------------------------------------------------------------------------------------------------

% For +z face we have Jx,Jy, Mx,My according to J = z x H ; and M = -z x E;
fdjxzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current x component array calculated at +z surface
fdjyzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current y component array calculated at +z surface
fdmxzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current x component array calculated at +z surface
fdmyzp=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current y component array calculated at +z surface

% For -z face we have Jx,Jy, Mx,My according to J = -z x H ; and M = -(-z) x E;
fdjxzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current x component array calculated at -z surface
fdjyzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current y component array calculated at -z surface
fdmxzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current x component array calculated at -z surface
fdmyzn=zeros(size(freq,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current y component array calculated at -z surface

%--------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------End of frequency domain electric and magnetic current arrays--------------------------------

%-------------------------------------End of Near Field Far Field - Frequency Domain Transformation Arrays-----------------------------




tic
%%
for n=1:T
    
    Hxincident_past = Hxincident_current; Hyincident_past = Hyincident_current; Hzincident_past = Hzincident_current;
    Exincident_past = Exincident_current; Eyincident_past = Eyincident_current; Ezincident_past = Ezincident_current;
    
    Exincident_current = Exincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEx).*exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    Eyincident_current = Eyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEy).*exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    Ezincident_current = Ezincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEz).*exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    Hxincident_current = Hxincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHx).*exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    Hyincident_current = Hyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHy).*exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    Hzincident_current = Hzincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHz).*exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);
    
    
    
    Hx(:,:,:) = C1Hx.*Hx+C2Hx.*(Ey(1:nxtot+1,1:nytot,2:nztot+1)-Ey(1:nxtot+1,1:nytot,1:nztot)) ...
        + C3Hx.*(Ez(1:nxtot+1,2:nytot+1,1:nztot)-Ez(1:nxtot+1,1:nytot,1:nztot));
    
    Hy(:,:,:) = C1Hy.*Hy+C2Hy.*(Ez(2:nxtot+1,1:nytot+1,1:nztot)-Ez(1:nxtot,1:nytot+1,1:nztot)) ...
        + C3Hy.*(Ex(1:nxtot,1:nytot+1,2:nztot+1)-Ex(1:nxtot,1:nytot+1,1:nztot));
    
    Hz(:,:,:) = C1Hz.*Hz+C2Hz.*(Ex(1:nxtot,2:nytot+1,1:nztot+1)-Ex(1:nxtot,1:nytot,1:nztot+1))  ...
        + C3Hz.*(Ey(2:nxtot+1,1:nytot,1:nztot+1)-Ey(1:nxtot,1:nytot,1:nztot+1));
    
    Hx = Hx + CHx_current.*Hxincident_current + CHx_past.*Hxincident_past;
    Hy = Hy + CHy_current.*Hyincident_current + CHy_past.*Hyincident_past;
    Hz = Hz + CHz_current.*Hzincident_current + CHz_past.*Hzincident_past;
    % Calculation of Auxiliary Fields and Application of CPML to Magnetic Field Components
    %-------------------------------------------------------------------xleft-----------------------------------------------------------------------------------
    for ii = 1: pml_thickness
        Psi_hy_ez_xl(ii,:,:)=bxml(ii)*Psi_hy_ez_xl(ii,:,:)+axml(ii)*(Ez(ii+1,:,:)-Ez(ii,:,:));  % Recursive convolution
        Psi_hz_ey_xl(ii,:,:)=bxml(ii)*Psi_hz_ey_xl(ii,:,:)+axml(ii)*(Ey(ii+1,:,:)-Ey(ii,:,:));  % Recursive convolution
    end
    
    Hy(1:pml_thickness,:,:) = Hy(1:pml_thickness,:,:)+CPsi_hy_ez_xl(:,:,:) .* Psi_hy_ez_xl(:,:,:);
    Hz(1:pml_thickness,:,:) = Hz(1:pml_thickness,:,:)+CPsi_hz_ey_xl(:,:,:) .* Psi_hz_ey_xl(:,:,:);
    
    %-------------------------------------------------------------------xright-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_hy_ez_xr(ii,:,:) = bxmr(ii)*Psi_hy_ez_xr(ii,:,:)+ axmr(ii)*(Ez(ii+nxtot-pml_thickness+1,:,:)-Ez(ii+nxtot-pml_thickness,:,:)); % Recursive convolution
        Psi_hz_ey_xr(ii,:,:) = bxmr(ii)*Psi_hz_ey_xr(ii,:,:)+ axmr(ii)*(Ey(ii+nxtot-pml_thickness+1,:,:)-Ey(ii+nxtot-pml_thickness,:,:)); % Recursive convolution
    end
    
    Hy(nxtot-pml_thickness+1:nxtot,:,:) = Hy(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_hy_ez_xr(:,:,:).*Psi_hy_ez_xr(:,:,:);
    Hz(nxtot-pml_thickness+1:nxtot,:,:) = Hz(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_hz_ey_xr(:,:,:).*Psi_hz_ey_xr(:,:,:);
    
    %-------------------------------------------------------------------yleft-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_hz_ex_yl(:,jj,:) = byml(jj) * Psi_hz_ex_yl(:,jj,:)+ayml(jj)*(Ex(:,jj+1,:)-Ex(:,jj,:)); % Recursive convolution
        Psi_hx_ez_yl(:,jj,:) = byml(jj) * Psi_hx_ez_yl(:,jj,:)+ayml(jj)*(Ez(:,jj+1,:)-Ez(:,jj,:)); % Recursive convolution
    end
    Hz(:,1:pml_thickness,:)= Hz(:,1:pml_thickness,:)+CPsi_hz_ex_yl(:,:,:).*Psi_hz_ex_yl(:,:,:);
    Hx(:,1:pml_thickness,:)= Hx(:,1:pml_thickness,:)+CPsi_hx_ez_yl(:,:,:).*Psi_hx_ez_yl(:,:,:);
    %-------------------------------------------------------------------yright----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        
        Psi_hz_ex_yr(:,jj,:) = bymr(jj)*Psi_hz_ex_yr(:,jj,:)+aymr(jj)*(Ex(:,jj+nytot-pml_thickness+1,:)-Ex(:,jj+nytot-pml_thickness,:)); % Recursive convolution
        Psi_hx_ez_yr(:,jj,:) = bymr(jj)*Psi_hx_ez_yr(:,jj,:) + aymr(jj)*(Ez(:,jj+nytot-pml_thickness+1,:)-Ez(:,jj+nytot-pml_thickness,:)); % Recursive convolution
    end
    Hz(:,nytot-pml_thickness+1:nytot,:)=Hz(:,nytot-pml_thickness+1:nytot,:)+CPsi_hz_ex_yr(:,:,:).*Psi_hz_ex_yr(:,:,:);
    Hx(:,nytot-pml_thickness+1:nytot,:)=Hx(:,nytot-pml_thickness+1:nytot,:)+CPsi_hx_ez_yr(:,:,:).*Psi_hx_ez_yr(:,:,:);
    %-------------------------------------------------------------------zleft----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_hx_ey_zl(:,:,kk) = bzml(kk)*Psi_hx_ey_zl(:,:,kk)+azml(kk)*(Ey(:,:,kk+1)-Ey(:,:,kk)); % Recursive convolution
        Psi_hy_ex_zl(:,:,kk) = bzml(kk)*Psi_hy_ex_zl(:,:,kk)+azml(kk)*(Ex(:,:,kk+1)-Ex(:,:,kk)); % Recursive convolution
    end
    
    Hx(:,:,1:pml_thickness) = Hx(:,:,1:pml_thickness)+CPsi_hx_ey_zl(:,:,:).*Psi_hx_ey_zl(:,:,:); % Recursive convolution
    Hy(:,:,1:pml_thickness) = Hy(:,:,1:pml_thickness)+CPsi_hy_ex_zl(:,:,:).*Psi_hy_ex_zl(:,:,:);    % Recursive convolution
    %-------------------------------------------------------------------zright----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_hx_ey_zr(:,:,kk) = bzmr(kk) * Psi_hx_ey_zr(:,:,kk)+azmr(kk)*(Ey(:,:,nztot-pml_thickness+1+kk)-Ey(:,:,nztot-pml_thickness+kk)); % Recursive convolution
        Psi_hy_ex_zr(:,:,kk) = bzmr(kk) * Psi_hy_ex_zr(:,:,kk)+azmr(kk)*(Ex(:,:,nztot-pml_thickness+1+kk)-Ex(:,:,nztot-pml_thickness+kk)); % Recursive convolution
    end
    
    Hx(:,:,nztot-pml_thickness+1:nztot) = Hx(:,:,nztot-pml_thickness+1:nztot)+CPsi_hx_ey_zr(:,:,:).*Psi_hx_ey_zr(:,:,:); % Recursive convolution
    Hy(:,:,nztot-pml_thickness+1:nztot) = Hy(:,:,nztot-pml_thickness+1:nztot)+CPsi_hy_ex_zr(:,:,:).*Psi_hy_ex_zr(:,:,:); % Recursive convolution
    %-------------------------------------------------------------------end------------------------------------------------------------------------------------
    
    
    Ex(1:nxtot,2:nytot,2:nztot) = C1Ex(1:nxtot,2:nytot,2:nztot).*Ex(1:nxtot,2:nytot,2:nztot) ...
        + C2Ex(1:nxtot,2:nytot,2:nztot).*(Hz(1:nxtot,2:nytot,2:nztot)-Hz(1:nxtot,1:nytot-1,2:nztot)) ...
        + C3Ex(1:nxtot,2:nytot,2:nztot).*(Hy(1:nxtot,2:nytot,2:nztot)-Hy(1:nxtot,2:nytot,1:nztot-1));
    
    Ey(2:nxtot,1:nytot,2:nztot)=C1Ey(2:nxtot,1:nytot,2:nztot).*Ey(2:nxtot,1:nytot,2:nztot) ...
        + C2Ey(2:nxtot,1:nytot,2:nztot).*(Hx(2:nxtot,1:nytot,2:nztot)-Hx(2:nxtot,1:nytot,1:nztot-1)) ...
        + C3Ey(2:nxtot,1:nytot,2:nztot).*(Hz(2:nxtot,1:nytot,2:nztot)-Hz(1:nxtot-1,1:nytot,2:nztot));
    
    Ez(2:nxtot,2:nytot,1:nztot)=C1Ez(2:nxtot,2:nytot,1:nztot).*Ez(2:nxtot,2:nytot,1:nztot) ...
        + C2Ez(2:nxtot,2:nytot,1:nztot).*(Hy(2:nxtot,2:nytot,1:nztot)-Hy(1:nxtot-1,2:nytot,1:nztot)) ...
        + C3Ez(2:nxtot,2:nytot,1:nztot).*(Hx(2:nxtot,2:nytot,1:nztot)-Hx(2:nxtot,1:nytot-1,1:nztot));
    
    Ex = Ex + CEx_current.* Exincident_current + CEx_past.* Exincident_past;
    Ey = Ey + CEy_current.* Eyincident_current + CEy_past.* Eyincident_past;
    Ez = Ez + CEz_current.* Ezincident_current + CEz_past .* Ezincident_past;
    
    %Calculation of Auxiliary Fields and Application of CPML to Electric Field Components
    %-------------------------------------------------------------------xleft-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_ey_hz_xl(ii,:,:) = bxl(ii)*Psi_ey_hz_xl(ii,:,:)+axl(ii)*(Hz(ii+1,:,:)-Hz(ii,:,:)); % Recursive convolution
        Psi_ez_hy_xl(ii,:,:) = bxl(ii)*Psi_ez_hy_xl(ii,:,:)+axl(ii)*(Hy(ii+1,:,:)-Hy(ii,:,:)); % Recursive convolution
    end
    Ey(2:pml_thickness+1,:,:)=Ey(2:pml_thickness+1,:,:)+CPsi_ey_hz_xl.*Psi_ey_hz_xl;
    Ez(2:pml_thickness+1,:,:)=Ez(2:pml_thickness+1,:,:)+CPsi_ez_hy_xl.*Psi_ez_hy_xl;
    %-------------------------------------------------------------------xright-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_ey_hz_xr(ii,:,:) = bxr(ii)*Psi_ey_hz_xr(ii,:,:)+axr(ii)*(Hz(ii+nxtot-pml_thickness,:,:)-Hz(ii+nxtot-pml_thickness-1,:,:)); % Recursive convolution
        Psi_ez_hy_xr(ii,:,:) = bxr(ii)*Psi_ez_hy_xr(ii,:,:)+axr(ii)*(Hy(ii+nxtot-pml_thickness,:,:)-Hy(ii+nxtot-pml_thickness-1,:,:)); % Recursive convolution
    end
    Ey(nxtot-pml_thickness+1:nxtot,:,:) = Ey(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_ey_hz_xr.*Psi_ey_hz_xr;
    Ez(nxtot-pml_thickness+1:nxtot,:,:) = Ez(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_ez_hy_xr.*Psi_ez_hy_xr;
    %-------------------------------------------------------------------yleft-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_ez_hx_yl(:,jj,:) = byl(jj) * Psi_ez_hx_yl(:,jj,:)+ayl(jj)*(Hx(:,jj+1,:)-Hx(:,jj,:)); % Recursive convolution
        Psi_ex_hz_yl(:,jj,:) = byl(jj) * Psi_ex_hz_yl(:,jj,:)+ayl(jj)*(Hz(:,jj+1,:)-Hz(:,jj,:)); % Recursive convolution
    end
    Ez(:,2:pml_thickness+1,:) = Ez(:,2:pml_thickness+1,:)+CPsi_ez_hx_yl.*Psi_ez_hx_yl;
    Ex(:,2:pml_thickness+1,:) = Ex(:,2:pml_thickness+1,:)+CPsi_ex_hz_yl.*Psi_ex_hz_yl;
    %-------------------------------------------------------------------yright-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_ez_hx_yr(:,jj,:) = byr(jj)*Psi_ez_hx_yr(:,jj,:)+ayr(jj)*(Hx(:,jj+nytot-pml_thickness,:)-Hx(:,jj+nytot-pml_thickness-1,:)); % Recursive convolution
        Psi_ex_hz_yr(:,jj,:) = byr(jj)*Psi_ex_hz_yr(:,jj,:)+ayr(jj)*(Hz(:,jj+nytot-pml_thickness,:)-Hz(:,jj+nytot-pml_thickness-1,:)); % Recursive convolution
    end
    Ez(:,nytot-pml_thickness+1:nytot,:) = Ez(:,nytot-pml_thickness+1:nytot,:)+CPsi_ez_hx_yr.*Psi_ez_hx_yr;
    Ex(:,nytot-pml_thickness+1:nytot,:) = Ex(:,nytot-pml_thickness+1:nytot,:)+CPsi_ex_hz_yr.*Psi_ex_hz_yr;
    %-------------------------------------------------------------------zleft-----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_ex_hy_zl(:,:,kk) = bzl(kk)*Psi_ex_hy_zl(:,:,kk)+azl(kk)*(Hy(:,:,kk+1)-Hy(:,:,kk)); % Recursive convolution
        Psi_ey_hx_zl(:,:,kk) = bzl(kk)*Psi_ey_hx_zl(:,:,kk)+azl(kk)*(Hx(:,:,kk+1)-Hx(:,:,kk)); % Recursive convolution
    end
    
    Ex(:,:,2:pml_thickness+1) = Ex(:,:,2:pml_thickness+1)+CPsi_ex_hy_zl.*Psi_ex_hy_zl;
    Ey(:,:,2:pml_thickness+1) = Ey(:,:,2:pml_thickness+1)+CPsi_ey_hx_zl.*Psi_ey_hx_zl;
    %-------------------------------------------------------------------zright-----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_ex_hy_zr(:,:,kk) = bzr(kk)*Psi_ex_hy_zr(:,:,kk)+azr(kk)*(Hy(:,:,kk+nztot-pml_thickness)-Hy(:,:,kk+nztot-pml_thickness-1)); % Recursive convolution
        Psi_ey_hx_zr(:,:,kk) = bzr(kk)*Psi_ey_hx_zr(:,:,kk)+azr(kk)*(Hx(:,:,kk+nztot-pml_thickness)-Hx(:,:,kk+nztot-pml_thickness-1)); % Recursive convolution
    end
    Ex(:,:,nztot-pml_thickness+1:nztot) = Ex(:,:,nztot-pml_thickness+1:nztot)+CPsi_ex_hy_zr.*Psi_ex_hy_zr;
    Ey(:,:,nztot-pml_thickness+1:nztot) = Ey(:,:,nztot-pml_thickness+1:nztot)+CPsi_ey_hx_zr.*Psi_ey_hx_zr;
    %-------------------------------------------------------------------end-----------------------------------------------------------------------------------
    %slice(Exincident_current,[1:32],[],[]);caxis([-Exincamp,Exincamp]);colorbar;title(['Iter number : ', num2str(n),' of ',num2str(T)]);view(45,45);drawnow;
    %subplot(311);
    %slice(Ezincident_current,[20],[20],[20]);caxis([-1,1]);colorbar;title(['Ex. Iter number : ', num2str(n),' of ',num2str(T)]);
    %subplot(312);
    %slice(Ey,[10:5:20],[],[]);caxis([-1,1]);shading interp;colorbar;title(['Ey. Iter number : ', num2str(n),' of ',num2str(T)]);
    %subplot(313);
    %slice(Ez,[10:5:20],[],[]);caxis([-1,1]);shading interp;
    %colorbar;title(['Ez. Iter number : ', num2str(n),' of ',num2str(T)]);
    %drawnow;%hold off;
    
    %% Calculate Magnetic and Electric Currents
    
    %------------------------- xp face magnetic and electric currents , we have J_y=-H_z , J_z=H_y , M_y=E_z , M_z=-E_y--------------------------------------
    
    tdmyxp(1,1,:,:) =  0.5*(Ez(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ez(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    tdmzxp(1,1,:,:) = -0.5*(Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    tdjyxp(1,1,:,:) =-0.25*(Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke) ...
        + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    
    tdjzxp(1,1,:,:) =0.25*(Hy(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hy(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1) ...
        + Hy (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hy (ffield_ie-1,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %------------------------ xn face magnetic and electric currents , we have J_y=H_z , J_z=-H_y , M_y=-E_z , M_z=+E_y-------------------------------------
    
    tdmyxn(1,1,:,:) = -0.5 * (Ez(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ez(ffield_is,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    tdmzxn(1,1,:,:) =  0.5 * (Ey(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ey(ffield_is,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    
    tdjyxn(1,1,:,:) = 0.25*(Hz(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hz(ffield_is,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke) ...
        + Hz (ffield_is-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_is-1,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    
    tdjzxn(1,1,:,:) =-0.25*(Hy(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hy(ffield_is,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1) ...
        + Hy (ffield_is-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hy (ffield_is-1,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %------------------------yp face magnetic and electric currents , we have J_x=H_z , J_z=-H_x , M_x=-E_z , M_z=E_x-------------------------------------
    
    tdmxyp(1,:,1,:) = -0.5*(Ez(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Ez(ffield_is+1:ffield_ie,ffield_je,ffield_ks:ffield_ke-1));
    tdmzyp(1,:,1,:) =  0.5*(Ex(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Ex(ffield_is:ffield_ie-1,ffield_je,ffield_ks+1:ffield_ke));
    
    tdjxyp(1,:,1,:) =  0.25*(Hz(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Hz(ffield_is:ffield_ie-1,ffield_je,ffield_ks+1:ffield_ke) ...
        + Hz (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks+1:ffield_ke));
    tdjzyp(1,:,1,:) =-0.25*(Hx(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Hx(ffield_is+1:ffield_ie,ffield_je,ffield_ks:ffield_ke-1) ...
        + Hx (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks:ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_je-1,ffield_ks:ffield_ke-1));
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %------------------------yn face magnetic and electric currents , we have J_x=-H_z , J_z=H_x , M_x=E_z , M_z=-E_x-------------------------------------
    %
    tdmxyn(1,:,1,:) =  0.5 * (Ez(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Ez(ffield_is+1:ffield_ie,ffield_js,ffield_ks:ffield_ke-1));
    tdmzyn(1,:,1,:) = -0.5 * (Ex(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Ex(ffield_is:ffield_ie-1,ffield_js,ffield_ks+1:ffield_ke));
    
    
    tdjzyn(1,:,1,:) = 0.25*(Hx(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Hx(ffield_is+1:ffield_ie,ffield_js,ffield_ks:ffield_ke-1) ...
        + Hx (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks:ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_js-1,ffield_ks:ffield_ke-1));
    
    tdjxyn(1,:,1,:) =-0.25*(Hz(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Hz(ffield_is:ffield_ie-1,ffield_js,ffield_ks+1:ffield_ke) ...
        + Hz (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks:ffield_ke-1) + Hz (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks+1:ffield_ke));
    
    %------------------------- zp face magnetic and electric currents , we have J_x=-H_y , J_y=H_x , M_x=E_y , M_y=-E_x--------------------------------------
    
    tdmxzp(1,:,:,1) =  0.5*(Ey(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Ey(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke));
    tdmyzp(1,:,:,1) = -0.5*(Ex(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Ex(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke));
    
    tdjyzp(1,:,:,1) = 0.25*(Hx(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke) ...
        + Hx (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke-1));
    
    tdjxzp(1,:,:,1) =-0.25*(Hy(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Hy(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke) ...
        + Hy (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke-1) + Hy (ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke-1));
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %------------------------- zn face magnetic and electric currents , we have J_x=H_y , J_y=-H_x , M_x=-E_y , M_y=E_x-------------------------------------
    
    tdmxzn(1,:,:,1) = -0.5 * (Ey(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Ey(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks));
    tdmyzn(1,:,:,1) =  0.5 * (Ex(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Ex(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks));
    
    tdjyzn(1,:,:,1) =-0.25*(Hx(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks) ...
        + Hx (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks-1)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks-1));
    
    tdjxzn(1,:,:,1) = 0.25*(Hy(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Hy(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks) ...
        + Hy (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks-1) + Hy (ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks-1));
    %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------Transform Magnetic and Electric currents to Frequency Domain------------------------------------------------------
    farfield_w=freq*2*pi;
    
    for ff=1:size(freq,2)
        kernelh = exp(-1i*farfield_w(ff)*(n-0.5)*delt)*delt; % Exponential term for DFT , here we take 0.5 time index diff. for electric currents since they
        % depend on magnetic fields
        
        % +- x face
        fdjyxp(ff,:,:,:) = fdjyxp(ff,:,:,:) + kernelh * tdjyxp(1,:,:,:);
        fdjzxp(ff,:,:,:) = fdjzxp(ff,:,:,:) + kernelh * tdjzxp(1,:,:,:);
        fdjzxn(ff,:,:,:) = fdjzxn(ff,:,:,:) + kernelh * tdjzxn(1,:,:,:);
        fdjyxn(ff,:,:,:) = fdjyxn(ff,:,:,:) + kernelh * tdjyxn(1,:,:,:);
        % +- y face
        fdjxyp(ff,:,:,:) = fdjxyp(ff,:,:,:) + kernelh * tdjxyp(1,:,:,:);
        fdjzyp(ff,:,:,:) = fdjzyp(ff,:,:,:) + kernelh * tdjzyp(1,:,:,:);
        fdjxyn(ff,:,:,:) = fdjxyn(ff,:,:,:) + kernelh * tdjxyn(1,:,:,:);
        fdjzyn(ff,:,:,:) = fdjzyn(ff,:,:,:) + kernelh * tdjzyn(1,:,:,:);
        % +- z face
        fdjyzp(ff,:,:,:) = fdjyzp(ff,:,:,:) + kernelh * tdjyzp(1,:,:,:);
        fdjxzp(ff,:,:,:) = fdjxzp(ff,:,:,:) + kernelh * tdjxzp(1,:,:,:);
        fdjyzn(ff,:,:,:) = fdjyzn(ff,:,:,:) + kernelh * tdjyzn(1,:,:,:);
        fdjxzn(ff,:,:,:) = fdjxzn(ff,:,:,:) + kernelh * tdjxzn(1,:,:,:);
        
        kernele = exp(-1i*farfield_w(ff)*n*delt)*delt;  % Exponential term for DFT , for Magnetic currents
        
        % +- x face
        fdmxyp(ff,:,:,:) = fdmxyp(ff,:,:,:) + kernele * tdmxyp(1,:,:,:);
        fdmzyp(ff,:,:,:) = fdmzyp(ff,:,:,:) + kernele * tdmzyp(1,:,:,:);
        fdmxyn(ff,:,:,:) = fdmxyn(ff,:,:,:) + kernele * tdmxyn(1,:,:,:);
        fdmzyn(ff,:,:,:) = fdmzyn(ff,:,:,:) + kernele * tdmzyn(1,:,:,:);
        % +- y face
        fdmyxp(ff,:,:,:) = fdmyxp(ff,:,:,:) + kernele * tdmyxp(1,:,:,:);
        fdmzxp(ff,:,:,:) = fdmzxp(ff,:,:,:) + kernele * tdmzxp(1,:,:,:);
        fdmyxn(ff,:,:,:) = fdmyxn(ff,:,:,:) + kernele * tdmyxn(1,:,:,:);
        fdmzxn(ff,:,:,:) = fdmzxn(ff,:,:,:) + kernele * tdmzxn(1,:,:,:);
        % +- z face
        fdmyzp(ff,:,:,:) = fdmyzp(ff,:,:,:) + kernele * tdmyzp(1,:,:,:);
        fdmxzp(ff,:,:,:) = fdmxzp(ff,:,:,:) + kernele * tdmxzp(1,:,:,:);
        fdmxzn(ff,:,:,:) = fdmxzn(ff,:,:,:) + kernele * tdmxzn(1,:,:,:);
        fdmyzn(ff,:,:,:) = fdmyzn(ff,:,:,:) + kernele * tdmyzn(1,:,:,:);
    end
    %-------------------------------------------------------------End of Transformation------------------------------------------------------
    display(['FDTD is Running at Time instant: ',num2str( n ),'/',num2str( T )]);
end
%---------------------------------------------------------------End of FDTD Loop----------------------------------------------------------
fprintf('\nCalculating Far Field Data...\n');
%% Radiated Power Calculation
%ffield.rpower=zeros(1,1);

% Radiated powers are calculated at each observation frequency from  0.5*real( conj(J)xM ), we split J and M into 6 face

% for ff=1:1;
%     power = 0;
%     power = delx*dely* sum(sum(sum(fdmyzp(ff,:,:,:).*conj(fdjxzp(ff,:,:,:)) - fdmxzp(ff,:,:,:).* conj(fdjyzp(ff,:,:,:)))));
%     power = power - delx*dely* sum(sum(sum(fdmyzn(ff,:,:,:) .* conj(fdjxzn(ff,:,:,:)) - fdmxzn(ff,:,:,:) .* conj(fdjyzn(ff,:,:,:)))));
%     power = power + delx*delz* sum(sum(sum(fdmxyp(ff,:,:,:) .* conj(fdjzyp(ff,:,:,:)) - fdmzyp(ff,:,:,:) .* conj(fdjxyp(ff,:,:,:)))));
%     power = power - delx*delz* sum(sum(sum(fdmxyn(ff,:,:,:) .* conj(fdjzyn(ff,:,:,:)) - fdmzyn(ff,:,:,:) .* conj(fdjxyn(ff,:,:,:)))));
%     power = power + dely*delz* sum(sum(sum(fdmzxp(ff,:,:,:) .* conj(fdjyxp(ff,:,:,:)) - fdmyxp(ff,:,:,:) .* conj(fdjzxp(ff,:,:,:)))));
%     power = power - dely*delz* sum(sum(sum(fdmzxn(ff,:,:,:) .* conj(fdjyxn(ff,:,:,:)) - fdmyxn(ff,:,:,:) .* conj(fdjzxn(ff,:,:,:)))));
%     ffield.rpower(ff) = 0.5 * real(power);
% end
power = zeros(1, size(freq,2));
w = 2 * pi * freq;
for n = 1:T
    power = power + -(sqrt(2*exp(1))/tau)*(n*delt - t0).*exp(-((n*delt - t0)/tau).^2) * exp(-1i*w*n*delt);
end
power = power* delt;
incidentwave.incident_power = (0.5/eta) * (incidentwave_ETheta^2 + incidentwave_EPhi^2) * abs(power).^2;

%% Far Field Calculation from Frequency Domain Electric and Magnetic Currents
ffield_plane='xy';
% xy plane theta=90 phi=(-180)-(180)
ffield_theta = zeros(ffield_nangles, 1);
ffield_phi   = zeros(ffield_nangles, 1);
ffield_theta=ffield_theta+inc_scat_angles(3)*pi/180;
if type_sim==2
    ffield_phi=pi/180*[-180:angle_step:179]'; % Phi (-180),(180)z
else
    ffield_phi=inc_scat_angles(4)*pi/180;% Phi (-180),(180)z
end

exp_jk_rpr = zeros(ffield_nangles,1);
delx_sinth_cosphi = zeros(ffield_nangles,1);
dely_sinth_sinphi = zeros(ffield_nangles,1);
delz_costh = zeros(ffield_nangles,1);
dely_delz_costh_sinphi = zeros(ffield_nangles,1);
dely_delz_sinth = zeros(ffield_nangles,1);
dely_delz_cosphi = zeros(ffield_nangles,1);
delx_delz_costh_cosphi = zeros(ffield_nangles,1);
delx_delz_sinth = zeros(ffield_nangles,1);
delx_delz_sinphi = zeros(ffield_nangles,1);
delx_dely_costh_cosphi = zeros(ffield_nangles,1);
delx_dely_costh_sinphi = zeros(ffield_nangles,1);
delx_dely_sinphi = zeros(ffield_nangles,1);
delx_dely_cosphi = zeros(ffield_nangles,1);

ffield_dirTheta = zeros(size(freq,2),ffield_nangles);
ffield_dir= zeros(size(freq,2),ffield_nangles);
ffield_dirPhi = zeros(size(freq,2),ffield_nangles);

delx_sinth_cosphi = delx*sin(ffield_theta).*cos(ffield_phi);
dely_sinth_sinphi = dely*sin(ffield_theta).*sin(ffield_phi);
delz_costh = delz*cos(ffield_theta);
dely_delz_costh_sinphi = dely*delz*cos(ffield_theta).*sin(ffield_phi);
dely_delz_sinth = dely*delz*sin(ffield_theta);
dely_delz_cosphi = dely*delz*cos(ffield_phi);
delx_delz_costh_cosphi = delx*delz*cos(ffield_theta).*cos(ffield_phi);
delx_delz_sinth = delx*delz*sin(ffield_theta);
delx_delz_sinphi = delx*delz*sin(ffield_phi);
delx_dely_costh_cosphi = delx*dely*cos(ffield_theta).*cos(ffield_phi);
delx_dely_costh_sinphi = delx*dely*cos(ffield_theta).*sin(ffield_phi);
delx_dely_sinphi = delx*dely*sin(ffield_phi);
delx_dely_cosphi = delx*dely*cos(ffield_phi);

ffield_cx = 0.5*(ffield_ie+ffield_is);
ffield_cy = 0.5*(ffield_je+ffield_js);
ffield_cz = 0.5*(ffield_ke+ffield_ks);

% calculate directivity
for ff=1:size(freq,2)
    
    k = 2*pi*freq(ff)*(mu*eps)^0.5;
    Ntheta = zeros(ffield_nangles,1);
    Ltheta = zeros(ffield_nangles,1);
    Nphi = zeros(ffield_nangles,1);
    Lphi = zeros(ffield_nangles,1);
    rpr = zeros(ffield_nangles,1);
    
    for yy = ffield_js:ffield_je-1
        for zz =ffield_ks:ffield_ke-1
            %----------------------------------------------------------------- +-x surfaces-----------------------------------------------------------------------------
            
            %----------------------------------------------------------------- +x surface-----------------------------------------------------------------------------
            R = (ffield_ie - ffield_cx)*delx_sinth_cosphi+(yy-ffield_cy+0.5)*dely_sinth_sinphi+(zz-ffield_cz+0.5)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta + (fdjyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdjzxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            Nphi = Nphi + (fdjyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdmzxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            
            Lphi = Lphi + (fdmyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
            %----------------------------------------------------------------- -x surface-----------------------------------------------------------------------------
            
            R = (ffield_is - ffield_cx)*delx_sinth_cosphi + (yy-ffield_cy+0.5)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta  + (fdjyxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdjzxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            Nphi = Nphi + (fdjyxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmyxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdmzxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            Lphi = Lphi + (fdmyxn(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
        end
    end
    %----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    %----------------------------------------------------------------- +-y surfaces-----------------------------------------------------------------------------
    for xx =ffield_is:ffield_ie-1
        for zz =ffield_ks:ffield_ke-1
            
            %------------------------------------------------------------------+y surface-----------------------------------------------------------------------------
            
            R = (xx - ffield_cx + 0.5)*delx_sinth_cosphi + (ffield_je-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta  + (fdjxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdjzyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdmzyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
            %------------------------------------------------------------------ -y surface-----------------------------------------------------------------------------
            R = (xx - ffield_cx + 0.5)*delx_sinth_cosphi + (ffield_js-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta + (fdjxyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdjzyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdmzyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxyn(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
        end
    end
    %----------------------------------------------------------------------------------------------------------------------------------------------------------
    for xx =ffield_is:ffield_ie-1
        for yy =ffield_js:ffield_je-1
            %------------------------------------------------------------------+-z surfaces-----------------------------------------------------------------------------
            
            %------------------------------------------------------------------ +z surface-----------------------------------------------------------------------------
            
            R = (xx-ffield_cx+0.5)*delx_sinth_cosphi + (yy - ffield_cy + 0.5)*dely_sinth_sinphi + (ffield_ke-ffield_cz)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta + (fdjxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdjyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi ...
                +fdjyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdmyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi) .*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi+ ...
                fdmyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
            
            %------------------------------------------------------------------ -z surface-----------------------------------------------------------------------------
            
            R = (xx-ffield_cx+0.5)*delx_sinth_cosphi + (yy - ffield_cy + 0.5)*dely_sinth_sinphi + (ffield_ks-ffield_cz)*delz_costh;
            
            exp_jk_rpr = exp(1i*k*R);
            
            Ntheta = Ntheta + (fdjxzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdjyzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi...
                +fdjyzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdmyzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi) .*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi...
                +fdmyzn(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
        end
    end
    ffield_dirTheta(ff,:)  = ...
        (k^2./(8*pi*eta*incidentwave.incident_power(ff))) ...
        .* (abs(Lphi+eta*Ntheta).^2);
    ffield_dirPhi(ff,:)    = ...
        (k^2./(8*pi*eta*incidentwave.incident_power(ff))) ...
        .* (abs(Ltheta-eta*Nphi).^2);
    
end


step_size=10;
Nrings=4;
line_style1='b-';
line_style2='r--';
scale_type='dB';
plot_type='D';

% xy plane

const_theta = 90; % used for plot


% Defaults for this blog post
width = 3.5;     % Width in inches
height = 3.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 9;      % Fontsize
lw = 2;      % LineWidth
msz = 8;       % MarkerSize


scrsz = get(0,'ScreenSize');
f = figure;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100],'Color',[1 1 1]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

for mi=1:size(freq,2)
    
    set(gcf,'Color',[1 1 1]);
    pat1 = ffield_dirTheta(mi,:).';
    pat2 = ffield_dirPhi(mi,:).';
    
    
    pat1 = 10*log10(pat1);
    pat2 = 10*log10(pat2);
    
    patt1=pat1;
    patt2=pat2;
    
    %pattern_1=pat1-min(pat1);
    pattern_2=pat2;
    
    
    maxFDTD=max(pat1);
    
    
    maxval=maxFDTD;
    
    
    pat1=pat1-maxval;
    
    
    pat1(pat1<-30)=-30;
    
    pattern_1=pat1+30;
    
    if type_sim==1
        theta=inc_scat_angles(3);
    elseif type_sim==2
        theta=0:deg2rad(angle_step):2*pi-pi/180;
    end
    
    if type_sim==1
        p1=plot(freq(mi),(patt1),'k*');
        hold on;
    elseif type_sim==2
        p1=polar(theta',(pattern_1),'k');
    end
    
    
    set(p1,'linewidth',2);
    h2=legend('FDTD','location','southeast');
    set(h2,'FontSize',9);
    
    text(18*cos(pi/4),50*sin(pi/4),['{\sigma}_{max} :',num2str((ceil(max(maxFDTD)*10)/10)),' dB'],'fontweight','bold');
    
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print('Figure8','-depsc2','-r300');
    
end
