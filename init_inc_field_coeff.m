function [CEx_current, CEx_past, CEy_current, CEy_past, CEz_current, CEz_past, ...
   CHx_current, CHx_past, CHy_current, CHy_past, CHz_current, CHz_past ] = init_inc_field_coeff(eps_x, eps_y, eps_z, ...
   eps, mu_x, mu_y, mu_z, mu, delt, sigma_x, sigma_y, sigma_z, sigmam_x, sigmam_y, sigmam_z)

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

end