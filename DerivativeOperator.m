function [Dxe,Dxxe,Dxh,Dye,Dyye,Dyh,Dze,Dzze,Dzh] = DerivativeOperator(dx,dy,dz,Nx,Ny,Nz,Nxy,Nxyz)


%% Defining the derivative operators
% Nx=3;Ny=3;Nz=3;Nxy = Nx*Ny;Nxyz=Nx*Ny*Nz;WX=ones(Nxyz,1);%WY=WX;WZ=WX;
% dx = 1; dy = 1; dz = 1;
v1 = -ones(Nx,1);
v2 = ones(Nx,1);
Dxe = (1/dx) .* spdiags([v1,v2],[0,1],Nx,Nx);
Dxe = kron(speye(Ny*Nz),sparse(Dxe));
clear v1 v2

v1 = ones(Nx,1); v2 = -2*ones(Nx,1); v3 = ones(Nx,1);
Dxxe = (1/dx^2) .* spdiags([v1,v2,v3],-1:1,Nx,Nx);
Dxxe = kron(speye(Ny*Nz),sparse(Dxxe));
clear v1 v2 v3

v1 = -ones(Nxy,1);
v2 = ones(Nxy,1);
Dye = (1/dy) .* spdiags([v1,v2], [0,Nx],Nxy,Nxy);
Dye = kron(speye(Nz),sparse(Dye));
%     for i = Nxy+1:Nxy:Nxyz
%         v2(i:i+Nx-1) = 0;
%         %v1(i-Nx:i-1) = 0;
%     end
%     Dye = (1/dy)*spdiags([v1,v2],[0,Nx],Nxyz,Nxyz);
clear v1 v2

v1 = ones(Nxyz,1);
v2 = -2*ones(Nxyz,1);
v3 = ones(Nxyz,1);
%     for i = Nxy+1:Nxy:Nxyz
%         v3(i:i+Nx-1) = 0;
%         v1(i-Nx:i-1) = 0;
%     end
Dyye = (1/dy)^2*spdiags([v1,v2,v3],[-Nx,0,Nx],Nxy,Nxy);
Dyye = kron(speye(Nz),sparse(Dyye));
clear v1 v2 v3;

v1 = -ones(Nxyz,1);
v2 = ones(Nxyz,1);
Dze = (1/dz)*spdiags([v1,v2],[0,Nxy],Nxyz,Nxyz);

clear v1 v2

v1 = ones(Nxyz,1);
v2 = -2*ones(Nxyz,1);
v3 = ones(Nxyz,1);
Dzze = (1/dz)^2*spdiags([v1,v2,v3],[-Nxy,0,Nxy],Nxyz,Nxyz);

clear v1 v2 v3

v1 = -ones(Nx,1);
v2 = ones(Nx,1);
Dxh = (1/dx) .* spdiags([v1,v2],[-1,0],Nx,Nx);
Dxh = kron(speye(Ny*Nz),sparse(Dxh));
clear v1 v2
% v1 = -ones(Nxyz,1);
% v2 = ones(Nxyz,1);
% Dxh = (1/dx) .* spdiags([v1,v2],-1:0,Nxyz,Nxyz);
% clear v1 v2
% index = 1;
% for i = Nx+1:Nx:Nxyz
%     Dxh(i,Nx*index) = 0;
%     index = index+1;
% end


v1 = -ones(Nxy,1);
v2 = ones(Nxy,1);
Dyh = (1/dy) .* spdiags([v1,v2], [-Nx,0],Nxy,Nxy);
Dyh = kron(speye(Nz),sparse(Dyh));
% v1 = -ones(Nxyz,1);
% v2 = ones(Nxyz,1);
% for i = Nxy-Nx:Nxy:Nxyz
%     v1(i+1:i+Nx) = 0;
% end
% Dyh = (1/dy)*spdiags([v1,v2],[-Nx,0],Nxyz,Nxyz);

clear v1 v2;

v1 = -ones(Nxyz,1);
v2 = ones(Nxyz,1);
Dzh = (1/dz) .* spdiags([v1,v2],[-Nxy 0],Nxyz,Nxyz);
clear v1 v2
% v1 = ones(Nxyz,1);
% v2 = -2*ones(Nxyz,1);
% v3 = ones(Nxyz,1);
% Dzzh = (1/dz)^2*spdiags([v1,v2,v3],[-Nxy,0,Nxy],Nxyz,Nxyz);

end