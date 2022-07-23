% Reading a free-form NASTRAN file and creating node and element list

% Inititalising MATLAB
function generateNodeListElementListFromNASTRAN_3D(file_name)
% clc; close all; clear variables;

fid = fopen(file_name);
% skipping first 4 lines
for i = 1:4
    str = fgetl(fid);
end
% Reading the node points
nodeList = [];
str = fgetl(fid);
while(~strcmp(str,'$ Element data section'))
    str = regexp(str, ',' , 'split');
    node = [str2double(str{4}), str2double(str{5}), str2double(str{6})];
    nodeList = [nodeList;node];
    str = fgetl(fid);
end

elementList = [];
str = fgetl(fid);
while(~strcmp(str,'$ Property data section'))
    str = regexp(str, ',' , 'split');
    element = [str2double(str{3}), str2double(str{4}), str2double(str{5}), ...
        str2double(str{6}), str2double(str{7}), str2double(str{8}), ...
        str2double(str{9}),str2double(str{10}), str2double(str{11})];
    elementList = [elementList; element];
    str = fgetl(fid);
end
fclose(fid);
save nodeList.txt nodeList -ascii
save elementList.txt elementList -ascii
end