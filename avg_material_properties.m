function [eps_x, mu_x, sigma_x, sigmam_x, eps_y, mu_y, sigma_y, sigmam_y,...
    eps_z, mu_z, sigma_z, sigmam_z] = avg_material_properties(materials, nxtot, nytot, nztot, sim_domian_materials,...
    eps_x, eps_y, eps_z,...
    mu_x, mu_y, mu_z, ...
    sigma_x, sigma_y, sigma_z, ...
    sigmam_x, sigmam_y, sigmam_z)
for ii = 1:size(materials,2)
    tempepsr(ii)   = materials(ii).epsr;
    tempmur(ii)    = materials(ii).mur;
    tempsigma(ii) = materials(ii).sigma;
    tempsigmam(ii) = materials(ii).sigmam;
end
tempmur(find(tempmur==0)) = 1e-20; % Avoid zero division
tempsigmam(find(tempsigmam==0)) = 1e-20; % Avoid zero division
%---------------------------------------------------------------------x components-----------------------------------------------------------------------
eps_x(1:nxtot,2:nytot,2:nztot)=0.25*(tempepsr(sim_domian_materials(1:nxtot,2:nytot,2:nztot))+tempepsr(sim_domian_materials(1:nxtot,1:nytot-1,2:nztot)) ...
    + tempepsr(sim_domian_materials(1:nxtot,2:nytot,1:nztot-1))+ tempepsr(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot-1)));

mu_x(2:nxtot,1:nytot,1:nztot)=2*(tempmur(sim_domian_materials(2:nxtot,1:nytot,1:nztot)).*tempmur(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot))) ...
    ./(tempmur(sim_domian_materials(2:nxtot,1:nytot,1:nztot))+tempmur(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot)));

sigma_x(1:nxtot,2:nytot,2:nztot)=0.25 * (tempsigma(sim_domian_materials(1:nxtot,2:nytot,2:nztot))+tempsigma(sim_domian_materials(1:nxtot,1:nytot-1,2:nztot)) ...
    + tempsigma(sim_domian_materials(1:nxtot,2:nytot,1:nztot-1))+tempsigma(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot-1)));

sigmam_x(2:nxtot,1:nytot,1:nztot)=2*(tempsigmam(sim_domian_materials(2:nxtot,1:nytot,1:nztot)).*tempsigmam(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot))) ...
    ./(tempsigmam(sim_domian_materials(2:nxtot,1:nytot,1:nztot))+tempsigmam(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot)));
%---------------------------------------------------------------------y components-----------------------------------------------------------------------
eps_y(2:nxtot,1:nytot,2:nztot)=0.25*(tempepsr(sim_domian_materials(2:nxtot,1:nytot,2:nztot))+tempepsr(sim_domian_materials(1:nxtot-1,1:nytot,2:nztot)) ...
    + tempepsr(sim_domian_materials(2:nxtot,1:nytot,1:nztot-1))+tempepsr(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot-1)));

mu_y(1:nxtot,2:nytot,1:nztot)=2*(tempmur(sim_domian_materials(1:nxtot,2:nytot,1:nztot)).*tempmur(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot))) ...
    ./(tempmur(sim_domian_materials(1:nxtot,2:nytot,1:nztot))+tempmur(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot)));

sigma_y(2:nxtot,1:nytot,2:nztot)=0.25*(tempsigma(sim_domian_materials(2:nxtot,1:nytot,2:nztot))+tempsigma(sim_domian_materials(1:nxtot-1,1:nytot,2:nztot)) ...
    + tempsigma(sim_domian_materials(2:nxtot,1:nytot,1:nztot-1))+tempsigma(sim_domian_materials(1:nxtot-1,1:nytot,1:nztot-1)));

sigmam_y(1:nxtot,2:nytot,1:nztot)=2*(tempsigmam(sim_domian_materials(1:nxtot,2:nytot,1:nztot)).*tempsigmam(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot))) ...
    ./(tempsigmam(sim_domian_materials(1:nxtot,2:nytot,1:nztot))+tempsigmam(sim_domian_materials(1:nxtot,1:nytot-1,1:nztot)));
%---------------------------------------------------------------------z components-----------------------------------------------------------------------
eps_z(2:nxtot,2:nytot,1:nztot)=0.25*(tempepsr(sim_domian_materials(2:nxtot,2:nytot,1:nztot))+tempepsr(sim_domian_materials(1:nxtot-1,2:nytot,1:nztot)) ...
    + tempepsr(sim_domian_materials(2:nxtot,1:nytot-1,1:nztot))+tempepsr(sim_domian_materials(1:nxtot-1,1:nytot-1,1:nztot)));

mu_z(1:nxtot,1:nytot,2:nztot)=2*(tempmur(sim_domian_materials(1:nxtot,1:nytot,2:nztot)).*tempmur(sim_domian_materials(1:nxtot,1:nytot,1:nztot-1))) ...
    ./(tempmur(sim_domian_materials(1:nxtot,1:nytot,2:nztot))+tempmur(sim_domian_materials(1:nxtot,1:nytot,1:nztot-1)));

sigma_z(2:nxtot,2:nytot,1:nztot)=0.25*(tempsigma(sim_domian_materials(2:nxtot,2:nytot,1:nztot))+tempsigma(sim_domian_materials(1:nxtot-1,2:nytot,1:nztot)) ...
    + tempsigma(sim_domian_materials(2:nxtot,1:nytot-1,1:nztot))+tempsigma(sim_domian_materials(1:nxtot-1,1:nytot-1,1:nztot)));

sigmam_z(1:nxtot,1:nytot,2:nztot)=2*(tempsigmam(sim_domian_materials(1:nxtot,1:nytot,2:nztot)).*tempsigmam(sim_domian_materials(1:nxtot,1:nytot,1:nztot-1))) ...
    ./(tempsigmam(sim_domian_materials(1:nxtot,1:nytot,2:nztot))+tempsigmam(sim_domian_materials(1:nxtot,1:nytot,1:nztot-1)));

figure('Color', 'w', 'units', 'normalized', 'OuterPosition',[0 0 1 1]);
slice(eps_x,30,[],[]); shading interp; axis equal; colorbar; view([90,0]); drawnow;  
end