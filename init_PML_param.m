function [Psi_ey_hz_xl, Psi_ez_hy_xl, Psi_hy_ez_xl, Psi_hz_ey_xl, ...
    CPsi_ey_hz_xl, CPsi_ez_hy_xl, CPsi_hy_ez_xl, CPsi_hz_ey_xl, ...
    Psi_ey_hz_xr, Psi_ez_hy_xr, Psi_hy_ez_xr, Psi_hz_ey_xr, ...
    CPsi_ey_hz_xr, CPsi_ez_hy_xr, CPsi_hy_ez_xr, CPsi_hz_ey_xr, ...
    Psi_ex_hz_yl, Psi_ez_hx_yl, Psi_hx_ez_yl, Psi_hz_ex_yl, ...
    CPsi_ex_hz_yl, CPsi_ez_hx_yl, CPsi_hx_ez_yl, CPsi_hz_ex_yl, ...
    Psi_ex_hz_yr, Psi_ez_hx_yr, Psi_hx_ez_yr, Psi_hz_ex_yr,...
    CPsi_ex_hz_yr, CPsi_ez_hx_yr, CPsi_hx_ez_yr, CPsi_hz_ex_yr, ...
    Psi_ex_hy_zl, Psi_ey_hx_zl, Psi_hx_ey_zl, Psi_hy_ex_zl, ...
    CPsi_ex_hy_zl, CPsi_ey_hx_zl, CPsi_hx_ey_zl, CPsi_hy_ex_zl,...
    Psi_ex_hy_zr, Psi_ey_hx_zr, Psi_hx_ey_zr, Psi_hy_ex_zr, ...
    CPsi_ex_hy_zr, CPsi_ey_hx_zr, CPsi_hx_ey_zr, CPsi_hy_ex_zr, ...
    bxl, bxml, axl, axml, bxr, bxmr, axr, axmr, ...
    byl, byml, ayl, ayml, byr, bymr, ayr, aymr, ...
    bzl, bzml, azl, azml, bzr, bzmr, azr, azmr,...
    C3Ey, C2Ez, C2Hy, C3Hz, ...
    C2Ex, C3Ez, C3Hx, C2Hz, ...
    C3Ex, C2Ey, C2Hx, C3Hy] = init_PML_param(C3Ey, C2Ez, C2Hy, C3Hz, C2Ex, C3Ez, C3Hx, C2Hz, C3Ex, C2Ey, C2Hx, C3Hy,...
    delx, dely, delz, delt, nxtot, nytot, nztot, pml_order, pml_thickness, pml_sigmafactor, pml_kappamax, pml_alphamin, ...
    pml_alphamax)

c = 2.99792458e8;   % Speed of Light
mu = 4.0*pi*1.0e-7; % Free space permeability
eps = 1.0/c^2/mu;   % Free space permittivity

sigmaopt=(pml_order+1)/(150*pi*delx); % Optimal sigma used for grading conductivity in x direction
sigmamax=pml_sigmafactor*sigmaopt; % Maximum sigma used for grading conductivity in x direction

% xleft face sigma(x) and sigmam(x), kappa(x) , kappam(x), alpha(x) and alpham(x) grading

% The sigma(x),kappa(x),alpha(x) are used to update field components Hy and Hz , sigmam(x),kappam(x) and alpham(x)
% are used to update field components Ez and Ey.

% The positions of Hy and Hz in x direction has an offset of (i+1/2) , on the other hand Ez and Ey has an offset of (i)


% x-left side
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

% x-right side
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

%y-left side
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

% y-right side
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

% z-left side
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

% z-right side
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

% end of CPML
end