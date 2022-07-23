% Reading a free-form NASTRAN file and creating node and element list

% Inititalising MATLAB
function generateNodeListElementListFromNASTRAN_2D(file_name)
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
    element = [str2double(str{4}), str2double(str{5}), str2double(str{6})];
    elementList = [elementList; element];
    str = fgetl(fid);
end
fclose(fid);
save nodeList_2D.txt nodeList -ascii
save elementList_2D.txt elementList -ascii
end