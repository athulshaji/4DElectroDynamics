function [eps_x, eps_y, eps_z,...
    mu_x, mu_y, mu_z,...
    sigma_x, sigma_y, sigma_z,...
    sigmam_x, sigmam_y, sigmam_z,...
    sim_domian_materials] = create_simulation_obj(  nxtot,...
    nytot,...
    nztot,...
    delx,...
    dely,...
    delz,...
    domain_lx,...
    domain_ly,...
    domain_lz,...
    modified_node_list)


% eps_x=ones(nxtot,nytot,nztot);   % Relative electrical permittivity array for x directiom
% eps_y=ones(nxtot,nytot,nztot);   % Relative electrical permittivity array for y directiom
% eps_z=ones(nxtot,nytot,nztot);   % Relative electrical permittivity array for z directiom
eps_x=ones(nxtot,nytot+1,nztot+1);   % Relative electrical permittivity array for x directiom
eps_y=ones(nxtot+1,nytot,nztot+1);   % Relative electrical permittivity array for y directiom
eps_z=ones(nxtot+1,nytot+1,nztot);   % Relative electrical permittivity array for z directiom

% mu_x=ones(nxtot,nytot,nztot);      % Relative magnetic permeability array for x direction
% mu_y=ones(nxtot,nytot,nztot);      % Relative magnetic permeability array for y direction
% mu_z=ones(nxtot,nytot,nztot);      % Relative magnetic permeability array for z direction
mu_x=ones(nxtot+1,nytot,nztot);      % Relative magnetic permeability array for x direction
mu_y=ones(nxtot,nytot+1,nztot);      % Relative magnetic permeability array for y direction
mu_z=ones(nxtot,nytot,nztot+1);      % Relative magnetic permeability array for z direction

% sigma_x=zeros(nxtot,nytot,nztot); % Electrical conductivity array for x direction
% sigma_y=zeros(nxtot,nytot,nztot); % Electrical conductivity array for y direction
% sigma_z=zeros(nxtot,nytot,nztot); % Electrical conductivity array for z direction
sigma_x=zeros(nxtot,nytot+1,nztot+1); % Electrical conductivity array for x direction
sigma_y=zeros(nxtot+1,nytot,nztot+1); % Electrical conductivity array for y direction
sigma_z=zeros(nxtot+1,nytot+1,nztot); % Electrical conductivity array for z direction

% sigmam_x=zeros(nxtot,nytot,nztot); % Magnetic conductivity array for x direction
% sigmam_y=zeros(nxtot,nytot,nztot); % Magnetic conductivity array for y direction
% sigmam_z=zeros(nxtot,nytot,nztot); % Magnetic conductivity array for z direction
sigmam_x=zeros(nxtot+1,nytot,nztot); % Magnetic conductivity array for x direction
sigmam_y=zeros(nxtot,nytot+1,nztot); % Magnetic conductivity array for y direction
sigmam_z=zeros(nxtot,nytot,nztot+1); % Magnetic conductivity array for z direction

sim_domian_materials = ones(nxtot, nytot, nztot); % Simulation Material Array , initialized with free space
% [rectprismlx, rectprismly, rectprismlz, rectprismux, rectprismuy, rectprismuz, sim_domian_materials] = create_simulation_obj(nxtot, nytot, nztot,...
%     delx, dely, delz, lx, ly, lz, ux, uy, uz, domain_lx, domain_ly, domain_lz)
% domain_centerx= zeros(nxtot,nytot,nztot); % Array for storing x direction cell center coordinates
% domain_centery = zeros(nxtot,nytot,nztot); % Array for storing y direction cell center coordinates
% domain_centerz = zeros(nxtot,nytot,nztot); % Array for storing z direction cell center coordinates
%
% for ii = 1:nxtot
%     domain_centerx(ii,:,:) = ...
%         (ii - 0.5) * delx + domain_lx; % x direction cell center coordinates which is used for object creation
% end
% for jj = 1:nytot
%     domain_centery(:,jj,:) = ...
%         (jj - 0.5) * dely + domain_ly; % y direction cell center coordinates which is used for object creation
% end
% for kk = 1:nztot
%     domain_centerz(:,:,kk) = ...
%         (kk - 0.5) * delz + domain_lz; % z direction cell center coordinates which is used for object creation
% end

% Assign Rectangular Prisms Materials
for i = 1:length(modified_node_list(:,1,1))
    ax = round((modified_node_list(i, 1) - domain_lx)/delx) + 1;
    ay = round((modified_node_list(i, 2) - domain_ly)/dely) + 1;
    az = round((modified_node_list(i, 3) - domain_lz)/delz) + 1;
%     ax = round((modified_node_list(i, 1) - domain_lx)/delx);
%     ay = round((modified_node_list(i, 2) - domain_ly)/dely);
%     az = round((modified_node_list(i, 3) - domain_lz)/delz);
%     
    
    sim_domian_materials(ax, ay, az)=2;
    
    eps_x(ax,ay,az) = modified_node_list(i,4);
    eps_y(ax,ay,az) = modified_node_list(i,4);
    eps_z(ax,ay,az) = modified_node_list(i,4);
    
    mu_x(ax,ay,az) = modified_node_list(i,5);
    mu_y(ax,ay,az) = modified_node_list(i,5);
    mu_z(ax,ay,az) = modified_node_list(i,5);
    
    sigma_x(ax,ay,az) = modified_node_list(i,6);
    sigma_y(ax,ay,az) = modified_node_list(i,6);
    sigma_z(ax,ay,az) = modified_node_list(i,6);
    
    sigmam_x(ax,ay,az) = modified_node_list(i,7);
    sigmam_y(ax,ay,az) = modified_node_list(i,7);
    sigmam_z(ax,ay,az) = modified_node_list(i,7);
end
% rectprismlx = round((lx - domain_lx)/delx) + 1;
% rectprismly = round((ly - domain_ly)/dely) + 1;
% rectprismlz = round((lz - domain_lz)/delz) + 1;
%
% rectprismux = round((ux - domain_lx)/delx)+1;
% rectprismuy = round((uy - domain_ly)/dely)+1;
% rectprismuz = round((uz - domain_lz)/delz)+1;
%
% sim_domian_materials(rectprismlx:rectprismux-1, rectprismly:rectprismuy-1, rectprismlz:rectprismuz-1)=2;

end