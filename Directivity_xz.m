function Directivity_xz (  T, ...
                        tau, t0, delt, ...
                        freqs_of_interest, ...
                        eta, ffield_nangles, ...
                        incidentwave_ETheta, ...
                        incidentwave_EPhi, ...
                        ffield_ie, ffield_je, ffield_ke, ...
                        ffield_is, ffield_js, ffield_ks, ...
                        delx, dely, delz, ...
                        fdjyxp, fdjzxp, fdmyxp, fdmzxp, ...
                        fdjyxn, fdjzxn, fdmyxn, fdmzxn, ...
                        fdjxyn, fdjzyn, fdmxyn, fdmzyn, ...
                        fdjxyp, fdjzyp, fdmxyp, fdmzyp, ...
                        fdjxzp, fdjyzp, fdmxzp, fdmyzp, ...
                        fdjxzn, fdjyzn, fdmxzn, fdmyzn)
fprintf('\nCalculating Far Field Data in the XY plane...\n');
w = 2 * pi * freqs_of_interest;
% Radiated power
for mi=1:length(freqs_of_interest)
   power=0.0; 
   power=delx*dely*sum(sum(sum(fdmyzp(mi,:,:,:).*conj(fdjxzp(mi,:,:,:))...
       -...
       fdmxzp(mi,:,:,:).*conj(fdjyzp(mi,:,:,:)))));
   power=power-delx*dely*sum(sum(sum(fdmyzn(mi,:,:,:).*conj(fdjxzn(mi,:,:,:))...
       -...
       fdmxzn(mi,:,:,:).*conj(fdjyzn(mi,:,:,:)))));
   power=power+delx*delz*sum(sum(sum(fdmxyp(mi,:,:,:).*conj(fdjzyp(mi,:,:,:))...
       -...
       fdmzyp(mi,:,:,:).*conj(fdjxyp(mi,:,:,:)))));
   power=power-delx*delz*sum(sum(sum(fdmxyn(mi,:,:,:).*conj(fdjzyn(mi,:,:,:))...
       -...
       fdmzyn(mi,:,:,:).*conj(fdjxyn(mi,:,:,:)))));
   power=power+dely*delz*sum(sum(sum(fdmzxp(mi,:,:,:).*conj(fdjyxp(mi,:,:,:))...
       -...
       fdmyxp(mi,:,:,:).*conj(fdjzxp(mi,:,:,:)))));
   power=power-dely*delz*sum(sum(sum(fdmzxn(mi,:,:,:).*conj(fdjyxn(mi,:,:,:))...
       -...
       fdmyxn(mi,:,:,:).*conj(fdjzxn(mi,:,:,:)))));
   
   radiated_power(mi)=0.5*real(power);
end
%% Far Field Calculation from Frequency Domain Electric and Magnetic Currents
c = 2.99792458e8;   % Speed of Light
mu = 4.0*pi*1.0e-7; % Free space permeability
eps = 1.0/c^2/mu;   % Free space permittivity
ffield_plane='xy';
% xy plane theta=90 phi=(-180)-(180)
ffield_theta = zeros(ffield_nangles, 1);
ffield_phi   = zeros(ffield_nangles, 1);
% ffield_theta=ffield_theta+90*pi/180;

% ffield_theta=ffield_theta+90*pi/180;
% ffield_phi=pi/180*(linspace(-180,180,ffield_nangles))'; % Phi (-180),(180)z

ffield_theta=ffield_theta+pi/180*(linspace(-180,180,ffield_nangles))';%ffield_theta+90*pi/180;
ffield_phi=ffield_phi+0*pi/180;%pi/180*(linspace(-180,180,ffield_nangles))'; % Phi (-180),(180)z

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

ffield_dirTheta = zeros(size(freqs_of_interest,2),ffield_nangles);
ffield_dir= zeros(size(freqs_of_interest,2),ffield_nangles);
ffield_dirPhi = zeros(size(freqs_of_interest,2),ffield_nangles);

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
for ff=1:size(freqs_of_interest,2)
    
    k = 2*pi*freqs_of_interest(ff)*(mu*eps)^0.5;
    Ntheta = zeros(ffield_nangles,1);
    Ltheta = zeros(ffield_nangles,1);
    Nphi = zeros(ffield_nangles,1);
    Lphi = zeros(ffield_nangles,1);
    rpr = zeros(ffield_nangles,1);
    
    for yy = ffield_js:ffield_je-1
        for zz =ffield_ks:ffield_ke-1
            %----------------------------------------------------------------- +-x surfaces-----------------------------------------------------------------------------
            
            %----------------------------------------------------------------- +x surface-----------------------------------------------------------------------------
            %R = (ffield_ie - ffield_cx)*delx_sinth_cosphi+(yy-ffield_cy+0.5)*dely_sinth_sinphi+(zz-ffield_cz+0.5)*delz_costh;
            R = (ffield_ie - ffield_cx)*delx_sinth_cosphi+(yy-ffield_cy+dely/2)*dely_sinth_sinphi+(zz-ffield_cz+delz/2)*delz_costh;
            exp_jk_rpr = exp(-1i*k*R);
            
            Ntheta = Ntheta + (fdjyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdjzxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            Nphi = Nphi + (fdjyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_costh_sinphi ...
                - fdmzxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_sinth).*exp_jk_rpr;
            
            Lphi = Lphi + (fdmyxp(ff,1,yy-ffield_js+1,zz-ffield_ks+1).*dely_delz_cosphi).*exp_jk_rpr;
            
            %----------------------------------------------------------------- -x surface-----------------------------------------------------------------------------
            
            %R = (ffield_is - ffield_cx)*delx_sinth_cosphi + (yy-ffield_cy+0.5)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            R = (ffield_is - ffield_cx)*delx_sinth_cosphi + (yy-ffield_cy+dely/2)*dely_sinth_sinphi + (zz-ffield_cz+delz/2)*delz_costh;
            exp_jk_rpr = exp(-1i*k*R);
            
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
            
            %R = (xx - ffield_cx + 0.5)*delx_sinth_cosphi + (ffield_je-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            R = (xx - ffield_cx + delx/2)*delx_sinth_cosphi + (ffield_je-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+delz/2)*delz_costh;
            
            exp_jk_rpr = exp(-1i*k*R);
            
            Ntheta = Ntheta  + (fdjxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdjzyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_costh_cosphi ...
                - fdmzyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinth).*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxyp(ff,xx-ffield_is+1,1,zz-ffield_ks+1).*delx_delz_sinphi).*exp_jk_rpr;
            
            %------------------------------------------------------------------ -y surface-----------------------------------------------------------------------------
            %R = (xx - ffield_cx + 0.5)*delx_sinth_cosphi + (ffield_js-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+0.5)*delz_costh;
            R = (xx - ffield_cx + delx/2)*delx_sinth_cosphi + (ffield_js-ffield_cy)*dely_sinth_sinphi + (zz-ffield_cz+delz/2)*delz_costh;
            exp_jk_rpr = exp(-1i*k*R);
            
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
            
            %R = (xx-ffield_cx+0.5)*delx_sinth_cosphi + (yy - ffield_cy + 0.5)*dely_sinth_sinphi + (ffield_ke-ffield_cz)*delz_costh;
            R = (xx-ffield_cx+delx/2)*delx_sinth_cosphi + (yy - ffield_cy + dely/2)*dely_sinth_sinphi + (ffield_ke-ffield_cz)*delz_costh;
            exp_jk_rpr = exp(-1i*k*R);
            
            Ntheta = Ntheta + (fdjxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdjyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi).*exp_jk_rpr;
            
            Nphi = Nphi + (-fdjxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi ...
                +fdjyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
            
            Ltheta = Ltheta + (fdmxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_cosphi ...
                + fdmyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_costh_sinphi) .*exp_jk_rpr;
            
            Lphi = Lphi + (-fdmxzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_sinphi+ ...
                fdmyzp(ff,xx-ffield_is+1,yy-ffield_js+1,1).*delx_dely_cosphi).*exp_jk_rpr;
            
            %------------------------------------------------------------------ -z surface-----------------------------------------------------------------------------
            
            %R = (xx-ffield_cx+0.5)*delx_sinth_cosphi + (yy - ffield_cy + 0.5)*dely_sinth_sinphi + (ffield_ks-ffield_cz)*delz_costh;
            R = (xx-ffield_cx+delx/2)*delx_sinth_cosphi + (yy - ffield_cy + dely/2)*dely_sinth_sinphi + (ffield_ks-ffield_cz)*delz_costh;
            exp_jk_rpr = exp(-1i*k*R);
            
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
    D_theta(ff,:)  = ...
        (k^2./(8*pi*eta*radiated_power(ff))) ...
        .* (abs(eta*Ntheta+Lphi).^2);
    D_phi(ff,:)    = ...
        (k^2./(8*pi*eta*radiated_power(ff))) ...
        .* (abs(eta*Nphi-Ltheta).^2);
    
end
D=D_theta+D_phi;
phi=deg2rad(linspace(0,360,360));
figure,polarplot(phi,D);
disp(['Max. D_xz = ',num2str(max(D))]);
figure,polarplot(phi,10*log10(D));
rlim([min(10*log10(D)),max(10*log10(D))]);
disp(['Max. D_xz(dB) = ',num2str(10*log10(max(D)))]);
end