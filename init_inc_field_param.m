function [Exincamp, Eyincamp, Ezincamp, Hxincamp, Hyincamp, Hzincamp, ...
    d0, k_x, k_y, k_z, ...
    kr0, spatialshift, krEx, krEy, krEz, krHx, krHy, krHz] = init_inc_field_param(incidentwave_inctheta,...
    incidentwave_incphi, incidentwave_ETheta, incidentwave_EPhi, eta, nxtot, nytot, nztot, delx, dely, delz,...
    domain_lx, domain_ly, domain_lz, domain_ux, domain_uy, domain_uz)

c = 2.99792458e8;   % Speed of Light
inctheta = incidentwave_inctheta*pi/180;
incphi = incidentwave_incphi*pi/180;
ETheta = incidentwave_ETheta;
EPhi = incidentwave_EPhi;

Exincamp = ETheta * cos(inctheta) * cos(incphi)- EPhi * sin(incphi);
Eyincamp = ETheta * cos(inctheta) * sin(incphi)+ EPhi * cos(incphi);
Ezincamp = -ETheta * sin(inctheta);

% Exincamp = ETheta * cos(inctheta) * cos(incphi)- EPhi * sin(inctheta)*sin(incphi);
% Eyincamp = ETheta * cos(inctheta) * sin(incphi)+ EPhi * sin(inctheta)*cos(incphi);
% Ezincamp = -ETheta * sin(inctheta);

Hxincamp = (-1/eta)*(EPhi * cos(inctheta)*cos(incphi) + ETheta * sin(incphi));
Hyincamp = (-1/eta)*(EPhi * cos(inctheta)*sin(incphi) - ETheta * cos(incphi));
Hzincamp = (1/eta)*(EPhi * sin(inctheta));

% Hxincamp = (-1/eta)*(ETheta * cos(inctheta) * cos(incphi)- EPhi * sin(inctheta)*sin(incphi));
% Hyincamp = (-1/eta)*(ETheta * cos(inctheta) * sin(incphi)+ EPhi * sin(inctheta)*cos(incphi));
% Hzincamp = (1/eta)*(ETheta * sin(inctheta));

cordx = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing x coordinates
cordy = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing y coordinates
cordz = zeros(nxtot+1,nytot+1,nztot+1); % Array for storing z coordinates
% cordx = zeros(nxtot,nytot,nztot); % Array for storing x coordinates
% cordy = zeros(nxtot,nytot,nztot); % Array for storing y coordinates
% cordz = zeros(nxtot,nytot,nztot); % Array for storing z coordinates
for ii= 1:nxtot+1
    cordx(ii,:,:) =  domain_lx + (ii - 1) * delx;
end
for jj = 1:nytot+1
    cordy(:,jj,:) = domain_ly + (jj - 1) * dely;
end
for kk = 1:nztot+1
    cordz(:,:,kk) = domain_lz+(kk - 1) * delz;
end


d0 =[domain_lx domain_ly domain_lz; 
    domain_lx domain_ly domain_uz; 
    domain_lx domain_uy domain_lz; 
    domain_lx domain_uy domain_uz;
    domain_ux domain_ly domain_lz; 
    domain_ux domain_ly domain_uz; 
    domain_ux domain_uy domain_lz; 
    domain_ux domain_uy domain_uz;];

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
end