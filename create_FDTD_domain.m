function [nxtot, nytot, nztot, domain_lx, domain_ly, domain_lz, domain_ux, domain_uy, domain_uz, X, Y, Z] = cretae_FDTD_domain(delx, dely, delz, pml_thickness, lx, ly, lz, ux, uy, uz, freespace_cells)
domain_lx = lx - delx * freespace_cells -delx*pml_thickness; % Simulation Domain lower x coordinate with freespace and PML
domain_ly = ly - dely * freespace_cells -dely*pml_thickness; % Simulation Domain lower y coordinate with freespace and PML
domain_lz = lz- delz * freespace_cells - delz*pml_thickness; % Simulation Domain lower z coordinate with freespace and PML

domain_ux = ux+ delx * freespace_cells + delx*pml_thickness; % Simulation Domain upper x coordinate with freespace and PML
domain_uy = uy+ dely * freespace_cells + dely*pml_thickness; % Simulation Domain upper y coordinate with freespace and PML
domain_uz = uz+ delz * freespace_cells + delz*pml_thickness; % Simulation Domain upper z coordinate with freespace and PML

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
end