function [ freq,incidentwave_ETheta,incidentwave_EPhi,incidentwave_inctheta,incidentwave_incphi ] = wave_parameters( file_name )
fid = fopen(file_name);
str = fgetl(fid);
param = regexp(str, ',', 'split');
freq = str2double(param{1});
incidentwave_ETheta = str2double(param{2});
incidentwave_EPhi = str2double(param{3});
incidentwave_inctheta = str2double(param{4});
incidentwave_incphi = str2double(param{5});
fclose(fid);
end

