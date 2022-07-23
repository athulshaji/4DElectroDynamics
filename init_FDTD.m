function [Ex, Ey, Ez, Hx, Hy, Hz, Exincident_current, Eyincident_current, Ezincident_current, ...
    Exincident_past, Eyincident_past, Ezincident_past, ...
    Hxincident_current, Hyincident_current, Hzincident_current, ...
    Hxincident_past, Hyincident_past, Hzincident_past] = init_FDTD(nxtot, nytot, nztot)

% Ex=zeros(nxtot,nytot,nztot);  % Electric field x component array
% Ey=zeros(nxtot,nytot,nztot);  % Electric field y component array
% Ez=zeros(nxtot,nytot,nztot);  % Electric field z component array
Ex=zeros(nxtot,nytot+1,nztot+1);  % Electric field x component array
Ey=zeros(nxtot+1,nytot,nztot+1);  % Electric field y component array
Ez=zeros(nxtot+1,nytot+1,nztot);  % Electric field z component array

% Hx=zeros(nxtot,nytot,nztot);    % Magnetic field x component array
% Hy=zeros(nxtot,nytot,nztot);    % Magnetic field y component array
% Hz=zeros(nxtot,nytot,nztot);    % Magnetic field z component array
Hx=zeros(nxtot+1,nytot,nztot);    % Magnetic field x component array
Hy=zeros(nxtot,nytot+1,nztot);    % Magnetic field y component array
Hz=zeros(nxtot,nytot,nztot+1);    % Magnetic field z component array

% eps_x=ones(nxtot,nytot+1,nztot+1);   % Relative electrical permittivity array for x directiom
% eps_y=ones(nxtot+1,nytot,nztot+1);   % Relative electrical permittivity array for y directiom
% eps_z=ones(nxtot+1,nytot+1,nztot);   % Relative electrical permittivity array for z directiom
% 
% mu_x=ones(nxtot+1,nytot,nztot);      % Relative magnetic permeability array for x direction
% mu_y=ones(nxtot,nytot+1,nztot);      % Relative magnetic permeability array for y direction
% mu_z=ones(nxtot,nytot,nztot+1);      % Relative magnetic permeability array for z direction
% 
% 
% sigma_x=zeros(nxtot,nytot+1,nztot+1); % Electrical conductivity array for x direction
% sigma_y=zeros(nxtot+1,nytot,nztot+1); % Electrical conductivity array for y direction
% sigma_z=zeros(nxtot+1,nytot+1,nztot); % Electrical conductivity array for z direction
% 
% sigmam_x=zeros(nxtot+1,nytot,nztot); % Magnetic conductivity array for x direction
% sigmam_y=zeros(nxtot,nytot+1,nztot); % Magnetic conductivity array for y direction
% sigmam_z=zeros(nxtot,nytot,nztot+1); % Magnetic conductivity array for z direction

% Exincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ex component's current value
% Eyincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ey component's current value
% Ezincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ez component's current value
Exincident_current = zeros(nxtot,nytot+1,nztot+1); % Array for storing Incident wave Ex component's current value
Eyincident_current = zeros(nxtot+1,nytot,nztot+1); % Array for storing Incident wave Ey component's current value
Ezincident_current = zeros(nxtot+1,nytot+1,nztot); % Array for storing Incident wave Ez component's current value

% Exincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ex component's past value
% Eyincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ey component's past value
% Ezincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Ez component's past value
Exincident_past = zeros(nxtot,nytot+1,nztot+1); % Array for storing Incident wave Ex component's past value
Eyincident_past = zeros(nxtot+1,nytot,nztot+1); % Array for storing Incident wave Ey component's past value
Ezincident_past = zeros(nxtot+1,nytot+1,nztot); % Array for storing Incident wave Ez component's past value

% Hxincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hx component's current value
% Hyincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hy component's current value
% Hzincident_current = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hz component's current value
Hxincident_current = zeros(nxtot+1,nytot,nztot); % Array for storing Incident wave Hx component's current value
Hyincident_current = zeros(nxtot,nytot+1,nztot); % Array for storing Incident wave Hy component's current value
Hzincident_current = zeros(nxtot,nytot,nztot+1); % Array for storing Incident wave Hz component's current value

% Hxincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hx component's past value
% Hyincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hy component's past value
% Hzincident_past = zeros(nxtot,nytot,nztot); % Array for storing Incident wave Hz component's past value
Hxincident_past = zeros(nxtot+1,nytot,nztot); % Array for storing Incident wave Hx component's past value
Hyincident_past = zeros(nxtot,nytot+1,nztot); % Array for storing Incident wave Hy component's past value
Hzincident_past = zeros(nxtot,nytot,nztot+1); % Array for storing Incident wave Hz component's past value
end