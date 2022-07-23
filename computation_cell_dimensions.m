function [lx, ly, lz, ux, uy, uz, freespace_cells] = computation_cell_dimensions(file_name)
fid = fopen(file_name);
str = fgetl(fid);
d = regexp(str, ',', 'split');
lx = str2double(d{1}); % Rectangular Prism lower x coordinate
ly = str2double(d{2}); % Rectangular Prism lower y coordinate
lz = str2double(d{3}); % Rectangular Prism lower z coordinate

ux = str2double(d{4});  % Rectangular Prism upper x coordinate
uy = str2double(d{5});  % Rectangular Prism lower y coordinate
uz = str2double(d{6});  % Rectangular Prism lower z coordinate

str = fgetl(fid);
freespace_cells = str2double(str);
end