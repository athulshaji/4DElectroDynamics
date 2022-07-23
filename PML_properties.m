function [ pml_thickness, pml_order, pml_rdesired, pml_sigmafactor, pml_kappamax, pml_alphamin, pml_alphamax] = PML_properties( file_name )
fid = fopen(file_name);
str = fgetl(fid);
param = regexp(str, ',', 'split');
pml_thickness = str2double(param{1}); %10; % PML thickness
pml_order = str2double(param{2}); %3; % PML order used for polynomial grading;
pml_rdesired = exp(-16); % PML desired maximum reflection error
pml_sigmafactor = str2double(param{4}); %1.3; % PML sigma factor used for grading sigma
pml_kappamax = str2double(param{5}); %7;  % PML maximum kappa used for grading kappa
pml_alphamin = str2double(param{6}); %0; %PML minimum alpha used for grading alpha
pml_alphamax = str2double(param{7}); %0.05; %PML maximum alpha used for grading alpha
fclose(fid);
end

