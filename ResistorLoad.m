function [C1EzResistor,C2EzResistor,C3EzResistor] = ResistorLoad(R,eps_z,eps, delt,sigma_z,delx,dely,delz,portx,porty,portz)
CONST = delt*delz./(R*delx*dely);
C1EzResistor=(2*eps_z(portx,porty,portz)*eps-delt*sigma_z(portx,porty,portz)-CONST)./(2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST); % First coefficent used to update Ez
C2EzResistor=(2*delt)./((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*delx); % Second coefficent used to update Ez
C3EzResistor=-(2*delt)./((2*eps_z(portx,porty,portz)*eps+delt*sigma_z(portx,porty,portz)+CONST)*dely); % Third coefficent used to update Ez
end