function [freqs_of_interest,ffield_nangles, ffield_theta, ffield_phi, ...
    ffield_dirTheta, ffield_dir, ffield_dirPhi, ...
    ffield_is, ffield_js, ffield_ks, ...
    ffield_ie, ffield_je, ffield_ke, ...
    ffield_cx, ffield_cy, ffield_cz] = init_NF2FF(file_name, nxtot, nytot, nztot)

fid = fopen(file_name);
str = fgetl(fid);
param = str2double(regexp(str, ',', 'split'));
ffield_distfromboundary = param(1); % Huygen Surface's Distance from the PML Absorbing Boundary
angle_step = param(2);

ffield_nangles=round(360/angle_step);    % Far field number of angles


% ffield_plane='xy';  % Obsevation plane
ffield_theta=zeros(ffield_nangles,1); % Observation theta
ffield_phi=zeros(ffield_nangles,1); % Observation phi

str = fgetl(fid);
freqs_of_interest = str2double(regexp(str, ',', 'split'));
ffield_dirTheta = zeros(size(freqs_of_interest,2),ffield_nangles); % Theta directivity
ffield_dir = zeros(size(freqs_of_interest,2),ffield_nangles);
ffield_dirPhi = zeros(size(freqs_of_interest,2),ffield_nangles); % Phi directivity


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
end