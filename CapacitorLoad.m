function [C1EzCap,C2EzCap,C3EzCap] = CapacitorLoad(C,eps_z,eps, delt,sigma_z,delx,dely,delz,portx,porty,portz)
CONST = (2*C*delz)/(delx*dely);
C1EzCap=(2*eps_z(portx,porty,portz)*eps-delt*sigma_z(portx,porty,portz)+CONST)./(2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST); % First coefficent used to update Ez
C2EzCap=(2*delt)/((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*delx); % Second coefficent used to update Ez
C3EzCap=-(2*delt)/((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*dely); % Third coefficent used to update Ez
end