function [C1EzVS,C2EzVS,C3EzVS,C4EzVS] = VoltageSource(Rs,eps_z,eps, delt,sigma_z,delx,dely,delz,portx,porty,portz)
CONST = delt*delz/(Rs*delx*dely);
C1EzVS=(2*eps_z(portx,porty,portz)*eps-delt*sigma_z(portx,porty,portz)-CONST)./(2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST); % First coefficent used to update Ez
C2EzVS=(2*delt)./((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*delx); % Second coefficent used to update Ez
C3EzVS=-(2*delt)./((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*dely); % Third coefficent used to update Ez
C4EzVS=-(2*delt)./((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*(Rs*dely*delx));
end