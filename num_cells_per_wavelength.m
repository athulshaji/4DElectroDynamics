function N = num_cells_per_wavelength(file_name)
fid = fopen(file_name);
str = fgetl(fid);
N = str2double(str);
fclose(fid);
end