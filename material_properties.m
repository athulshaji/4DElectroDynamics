function [ materials ] = material_properties(file_name)
materials=[];   % Material Structure
fid = fopen(file_name);
num_of_materials = str2double(fgetl(fid));
for i = 1:num_of_materials
    str = fgetl(fid);
    m = str2double(regexp(str, ',', 'split'));
    materials(i).epsr = m(1); % Relative electrical permittivity
    materials(i).mur = m(2); % Relative magnetic permeability
    materials(i).sigma = m(3); % Material electric conductivity
    materials(i).sigmam = m(4); % Material magnetic conductivity   
end
fclose(fid);
end

