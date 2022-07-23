function modified_node_list = Preprocess(E, N, materials)
modified_node_list = zeros(length(N), 7);

for i = 1:length(E)
    for j = 2:9
        modified_node_list(E(i,j), :) = [N(E(i,j),:), materials(E(i,1)).epsr, materials(E(i,1)).mur,...
        materials(E(i,1)).sigma, materials(E(i,1)).sigmam];
    end
end
% for i = 1:length(N)
%     modified_node_list(i, :) = [N(i,:), material(1).epsr, material(1).mur,...
%         material(1).sigma, material(1).sigmam];
% end
% save modified_nodel_list.txt modified_node_list -ascii
end
% Load the node and element list
% E = load('data\elementList.txt');
% N = load('data\nodeList.txt');

% material = [];
% 
% material(1).epsr = 1e20;
% material(1).mur = 1;
% material(1).sigma = 0;
% material(1).sigmam = 0;
% 
% material(2).epsr = 1;
% material(2).mur = 1;
% material(2).sigma = 0;
% material(2).sigmar = 0;


% a = 51;
% startx = -1;
% starty = -1;
% startz = 0;
% 
% stopx = 1;
% stopy = 1;
% stopz = 1;
% 
% stepx = 0.1;
% stepy = 0.1;
% stepz = 0.1;
% 
% x = startx:stepx:stopx+stepx;
% y = starty:stepy:stopy+stepy;
% z = startz:stepz:stopz+stepz;
% a = length(x);
% b = length(y);
% c = length(z);
% epsr = ones(a, b, c);
% [X, Y, Z] = meshgrid(x, y, z);
% 
% for q = 1:length(N)
%     xq = modified_node_list(q,1);
%     yq = modified_node_list(q,2);
%     zq = modified_node_list(q,3);
%     nn = [x(1), y(1), z(1)];
%     old_dist = ((x(1) - xq)^2 + ...
%                     (y(1) - yq)^2 + ...
%                     (z(1) - zq)^2);
%     for i = 1:a
%         for j = 1:b
%             for k = 1:c
%                 dist = ((x(i) - xq)^2 + ...
%                     (y(j) - yq)^2 + ...
%                     (z(k) - zq)^2);
%                 if dist <= old_dist
%                     nn = [x(i), y(j), z(k)];
%                     ii = i;
%                     jj = j;
%                     kk = k;
%                     old_dist = dist;
%                 end
%             end
%         end
%     end   
%     epsr(ii, jj, kk) = modified_node_list(q, 4);
% end
% 
% figure,
% xslice = [];
% yslice = [];
% zslice = [0, 0.5, 1];
% slice(X, Y, Z, epsr, xslice, yslice, zslice);
% colorbar
% 
% figure, xslice = [0];
% yslice = [];
% zslice = [];
% slice(X, Y, Z, epsr, xslice, yslice, zslice);
% colorbar
% 
% figure, xslice = [];
% yslice = [0];
% zslice = [];
% slice(X, Y, Z, epsr, xslice, yslice, zslice);
% colorbar
