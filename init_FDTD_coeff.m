function [C1Ex, C2Ex, C3Ex, C1Ey, C2Ey, C3Ey, C1Ez, C2Ez, C3Ez, ...
    C1Hx, C2Hx, C3Hx, C1Hy, C2Hy, C3Hy, C1Hz, C2Hz, C3Hz] = init_FDTD_coeff(eps_x, eps_y, eps_z, ...
   eps, mu_x, mu_y, mu_z, mu, delt, delx, dely, delz, sigma_x, sigma_y, sigma_z, sigmam_x, sigmam_y, sigmam_z)
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
end