function FDTD(  T, ...
                nxtot,...
                nytot,...
                nztot,...
                t0, ...
                tau, ...
                delt, ...  
                pml_thickness, ...
                krEx, ...
                krEy, ...
                krEz, ...
                krHx, ...
                krHy, ...
                krHz, ...
                Exincamp, ...
                Eyincamp, ...
                Ezincamp, ...
                Hxincamp, ...
                Hyincamp, ...
                Hzincamp, ...
                Exincident_current, ...
                Eyincident_current, ...
                Ezincident_current, ...
                Hxincident_current, ...
                Hyincident_current, ...
                Hzincident_current, ...
                Exincident_past, ...
                Eyincident_past, ...
                Ezincident_past, ...
                Hxincident_past, ...
                Hyincident_past, ...
                Hzincident_past, ...
                CHx_current, CHy_current, CHz_current, ...
                CHx_past, CHy_past, CHz_past, ...
                CEx_current, CEy_current, CEz_current, ...
                CEx_past, CEy_past, CEz_past, ...
                Ex, Ey, Ez, Hx, Hy, Hz, ...
                C1Ex, C2Ex, C3Ex, C1Ey, C2Ey, C3Ey, C1Ez, C2Ez, C3Ez, ...
                C1Hx, C2Hx, C3Hx, C1Hy, C2Hy, C3Hy, C1Hz, C2Hz, C3Hz, ...
                Psi_ey_hz_xl, Psi_ez_hy_xl, Psi_hy_ez_xl, Psi_hz_ey_xl, ...
                CPsi_ey_hz_xl, CPsi_ez_hy_xl, CPsi_hy_ez_xl, CPsi_hz_ey_xl, ...
                Psi_ey_hz_xr, Psi_ez_hy_xr, Psi_hy_ez_xr, Psi_hz_ey_xr, ...
                CPsi_ey_hz_xr, CPsi_ez_hy_xr, CPsi_hy_ez_xr, CPsi_hz_ey_xr, ...
                Psi_ex_hz_yl, Psi_ez_hx_yl, Psi_hx_ez_yl, Psi_hz_ex_yl, ...
                CPsi_ex_hz_yl, CPsi_ez_hx_yl, CPsi_hx_ez_yl, CPsi_hz_ex_yl, ...
                Psi_ex_hz_yr, Psi_ez_hx_yr, Psi_hx_ez_yr, Psi_hz_ex_yr,...
                CPsi_ex_hz_yr, CPsi_ez_hx_yr, CPsi_hx_ez_yr, CPsi_hz_ex_yr, ...
                Psi_ex_hy_zl, Psi_ey_hx_zl, Psi_hx_ey_zl, Psi_hy_ex_zl, ...
                CPsi_ex_hy_zl, CPsi_ey_hx_zl, CPsi_hx_ey_zl, CPsi_hy_ex_zl,...
                Psi_ex_hy_zr, Psi_ey_hx_zr, Psi_hx_ey_zr, Psi_hy_ex_zr, ...
                CPsi_ex_hy_zr, CPsi_ey_hx_zr, CPsi_hx_ey_zr, CPsi_hy_ex_zr, ...
                bxl, bxml, axl, axml, bxr, bxmr, axr, axmr, ...
                byl, byml, ayl, ayml, byr, bymr, ayr, aymr, ...
                bzl, bzml, azl, azml, bzr, bzmr, azr, azmr)

fig = figure('Color', 'w', 'units', 'normalized', 'OuterPosition',[0 0 1 1]);
movieName = '3D_FDTD.avi';
vidObj = VideoWriter(movieName);
open(vidObj);
for n=1:T
    Hxincident_past = Hxincident_current;
    Hyincident_past = Hyincident_current; 
    Hzincident_past = Hzincident_current;
    Exincident_past = Exincident_current; 
    Eyincident_past = Eyincident_current; 
    Ezincident_past = Ezincident_current;
  
    Exincident_current = Exincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEx).*exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    Eyincident_current = Eyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEy).*exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    Ezincident_current = Ezincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEz).*exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    Hxincident_current = Hxincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHx).*exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    Hyincident_current = Hyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHy).*exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    Hzincident_current = Hzincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHz).*exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);
    
    
    
    Hx(:,:,:) = C1Hx.*Hx+C2Hx.*(Ey(1:nxtot+1,1:nytot,2:nztot+1)-Ey(1:nxtot+1,1:nytot,1:nztot)) ...
        + C3Hx.*(Ez(1:nxtot+1,2:nytot+1,1:nztot)-Ez(1:nxtot+1,1:nytot,1:nztot));
    
    Hy(:,:,:) = C1Hy.*Hy+C2Hy.*(Ez(2:nxtot+1,1:nytot+1,1:nztot)-Ez(1:nxtot,1:nytot+1,1:nztot)) ...
        + C3Hy.*(Ex(1:nxtot,1:nytot+1,2:nztot+1)-Ex(1:nxtot,1:nytot+1,1:nztot));
    
    Hz(:,:,:) = C1Hz.*Hz+C2Hz.*(Ex(1:nxtot,2:nytot+1,1:nztot+1)-Ex(1:nxtot,1:nytot,1:nztot+1))  ...
        + C3Hz.*(Ey(2:nxtot+1,1:nytot,1:nztot+1)-Ey(1:nxtot,1:nytot,1:nztot+1));
    
    Hx = Hx + CHx_current.*Hxincident_current + CHx_past.*Hxincident_past;
    Hy = Hy + CHy_current.*Hyincident_current + CHy_past.*Hyincident_past;
    Hz = Hz + CHz_current.*Hzincident_current + CHz_past.*Hzincident_past;
    % Calculation of Auxiliary Fields and Application of CPML to Magnetic Field Components
    %-------------------------------------------------------------------xleft-----------------------------------------------------------------------------------
    for ii = 1: pml_thickness
        Psi_hy_ez_xl(ii,:,:)=bxml(ii)*Psi_hy_ez_xl(ii,:,:)+axml(ii)*(Ez(ii+1,:,:)-Ez(ii,:,:));  % Recursive convolution
        Psi_hz_ey_xl(ii,:,:)=bxml(ii)*Psi_hz_ey_xl(ii,:,:)+axml(ii)*(Ey(ii+1,:,:)-Ey(ii,:,:));  % Recursive convolution
    end
    
    Hy(1:pml_thickness,:,:) = Hy(1:pml_thickness,:,:)+CPsi_hy_ez_xl(:,:,:) .* Psi_hy_ez_xl(:,:,:);
    Hz(1:pml_thickness,:,:) = Hz(1:pml_thickness,:,:)+CPsi_hz_ey_xl(:,:,:) .* Psi_hz_ey_xl(:,:,:);
    
    %-------------------------------------------------------------------xright-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_hy_ez_xr(ii,:,:) = bxmr(ii)*Psi_hy_ez_xr(ii,:,:)+ axmr(ii)*(Ez(ii+nxtot-pml_thickness+1,:,:)-Ez(ii+nxtot-pml_thickness,:,:)); % Recursive convolution
        Psi_hz_ey_xr(ii,:,:) = bxmr(ii)*Psi_hz_ey_xr(ii,:,:)+ axmr(ii)*(Ey(ii+nxtot-pml_thickness+1,:,:)-Ey(ii+nxtot-pml_thickness,:,:)); % Recursive convolution
    end
    
    Hy(nxtot-pml_thickness+1:nxtot,:,:) = Hy(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_hy_ez_xr(:,:,:).*Psi_hy_ez_xr(:,:,:);
    Hz(nxtot-pml_thickness+1:nxtot,:,:) = Hz(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_hz_ey_xr(:,:,:).*Psi_hz_ey_xr(:,:,:);
    
    %-------------------------------------------------------------------yleft-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_hz_ex_yl(:,jj,:) = byml(jj) * Psi_hz_ex_yl(:,jj,:)+ayml(jj)*(Ex(:,jj+1,:)-Ex(:,jj,:)); % Recursive convolution
        Psi_hx_ez_yl(:,jj,:) = byml(jj) * Psi_hx_ez_yl(:,jj,:)+ayml(jj)*(Ez(:,jj+1,:)-Ez(:,jj,:)); % Recursive convolution
    end
    Hz(:,1:pml_thickness,:)= Hz(:,1:pml_thickness,:)+CPsi_hz_ex_yl(:,:,:).*Psi_hz_ex_yl(:,:,:);
    Hx(:,1:pml_thickness,:)= Hx(:,1:pml_thickness,:)+CPsi_hx_ez_yl(:,:,:).*Psi_hx_ez_yl(:,:,:);
    %-------------------------------------------------------------------yright----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        
        Psi_hz_ex_yr(:,jj,:) = bymr(jj)*Psi_hz_ex_yr(:,jj,:)+aymr(jj)*(Ex(:,jj+nytot-pml_thickness+1,:)-Ex(:,jj+nytot-pml_thickness,:)); % Recursive convolution
        Psi_hx_ez_yr(:,jj,:) = bymr(jj)*Psi_hx_ez_yr(:,jj,:) + aymr(jj)*(Ez(:,jj+nytot-pml_thickness+1,:)-Ez(:,jj+nytot-pml_thickness,:)); % Recursive convolution
    end
    Hz(:,nytot-pml_thickness+1:nytot,:)=Hz(:,nytot-pml_thickness+1:nytot,:)+CPsi_hz_ex_yr(:,:,:).*Psi_hz_ex_yr(:,:,:);
    Hx(:,nytot-pml_thickness+1:nytot,:)=Hx(:,nytot-pml_thickness+1:nytot,:)+CPsi_hx_ez_yr(:,:,:).*Psi_hx_ez_yr(:,:,:);
    %-------------------------------------------------------------------zleft----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_hx_ey_zl(:,:,kk) = bzml(kk)*Psi_hx_ey_zl(:,:,kk)+azml(kk)*(Ey(:,:,kk+1)-Ey(:,:,kk)); % Recursive convolution
        Psi_hy_ex_zl(:,:,kk) = bzml(kk)*Psi_hy_ex_zl(:,:,kk)+azml(kk)*(Ex(:,:,kk+1)-Ex(:,:,kk)); % Recursive convolution
    end
    
    Hx(:,:,1:pml_thickness) = Hx(:,:,1:pml_thickness)+CPsi_hx_ey_zl(:,:,:).*Psi_hx_ey_zl(:,:,:); % Recursive convolution
    Hy(:,:,1:pml_thickness) = Hy(:,:,1:pml_thickness)+CPsi_hy_ex_zl(:,:,:).*Psi_hy_ex_zl(:,:,:);    % Recursive convolution
    %-------------------------------------------------------------------zright----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_hx_ey_zr(:,:,kk) = bzmr(kk) * Psi_hx_ey_zr(:,:,kk)+azmr(kk)*(Ey(:,:,nztot-pml_thickness+1+kk)-Ey(:,:,nztot-pml_thickness+kk)); % Recursive convolution
        Psi_hy_ex_zr(:,:,kk) = bzmr(kk) * Psi_hy_ex_zr(:,:,kk)+azmr(kk)*(Ex(:,:,nztot-pml_thickness+1+kk)-Ex(:,:,nztot-pml_thickness+kk)); % Recursive convolution
    end
    
    Hx(:,:,nztot-pml_thickness+1:nztot) = Hx(:,:,nztot-pml_thickness+1:nztot)+CPsi_hx_ey_zr(:,:,:).*Psi_hx_ey_zr(:,:,:); % Recursive convolution
    Hy(:,:,nztot-pml_thickness+1:nztot) = Hy(:,:,nztot-pml_thickness+1:nztot)+CPsi_hy_ex_zr(:,:,:).*Psi_hy_ex_zr(:,:,:); % Recursive convolution
    %-------------------------------------------------------------------end------------------------------------------------------------------------------------
    
    
    Ex(1:nxtot,2:nytot,2:nztot) = C1Ex(1:nxtot,2:nytot,2:nztot).*Ex(1:nxtot,2:nytot,2:nztot) ...
        + C2Ex(1:nxtot,2:nytot,2:nztot).*(Hz(1:nxtot,2:nytot,2:nztot)-Hz(1:nxtot,1:nytot-1,2:nztot)) ...
        + C3Ex(1:nxtot,2:nytot,2:nztot).*(Hy(1:nxtot,2:nytot,2:nztot)-Hy(1:nxtot,2:nytot,1:nztot-1));
    
    Ey(2:nxtot,1:nytot,2:nztot)=C1Ey(2:nxtot,1:nytot,2:nztot).*Ey(2:nxtot,1:nytot,2:nztot) ...
        + C2Ey(2:nxtot,1:nytot,2:nztot).*(Hx(2:nxtot,1:nytot,2:nztot)-Hx(2:nxtot,1:nytot,1:nztot-1)) ...
        + C3Ey(2:nxtot,1:nytot,2:nztot).*(Hz(2:nxtot,1:nytot,2:nztot)-Hz(1:nxtot-1,1:nytot,2:nztot));
    
    Ez(2:nxtot,2:nytot,1:nztot)=C1Ez(2:nxtot,2:nytot,1:nztot).*Ez(2:nxtot,2:nytot,1:nztot) ...
        + C2Ez(2:nxtot,2:nytot,1:nztot).*(Hy(2:nxtot,2:nytot,1:nztot)-Hy(1:nxtot-1,2:nytot,1:nztot)) ...
        + C3Ez(2:nxtot,2:nytot,1:nztot).*(Hx(2:nxtot,2:nytot,1:nztot)-Hx(2:nxtot,1:nytot-1,1:nztot));
    
    Ex = Ex + CEx_current.* Exincident_current + CEx_past.* Exincident_past;
    Ey = Ey + CEy_current.* Eyincident_current + CEy_past.* Eyincident_past;
    Ez = Ez + CEz_current.* Ezincident_current + CEz_past .* Ezincident_past;
    
    %Calculation of Auxiliary Fields and Application of CPML to Electric Field Components
    %-------------------------------------------------------------------xleft-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_ey_hz_xl(ii,:,:) = bxl(ii)*Psi_ey_hz_xl(ii,:,:)+axl(ii)*(Hz(ii+1,:,:)-Hz(ii,:,:)); % Recursive convolution
        Psi_ez_hy_xl(ii,:,:) = bxl(ii)*Psi_ez_hy_xl(ii,:,:)+axl(ii)*(Hy(ii+1,:,:)-Hy(ii,:,:)); % Recursive convolution
    end
    Ey(2:pml_thickness+1,:,:)=Ey(2:pml_thickness+1,:,:)+CPsi_ey_hz_xl.*Psi_ey_hz_xl;
    Ez(2:pml_thickness+1,:,:)=Ez(2:pml_thickness+1,:,:)+CPsi_ez_hy_xl.*Psi_ez_hy_xl;
    %-------------------------------------------------------------------xright-----------------------------------------------------------------------------------
    for ii = 1:pml_thickness
        Psi_ey_hz_xr(ii,:,:) = bxr(ii)*Psi_ey_hz_xr(ii,:,:)+axr(ii)*(Hz(ii+nxtot-pml_thickness,:,:)-Hz(ii+nxtot-pml_thickness-1,:,:)); % Recursive convolution
        Psi_ez_hy_xr(ii,:,:) = bxr(ii)*Psi_ez_hy_xr(ii,:,:)+axr(ii)*(Hy(ii+nxtot-pml_thickness,:,:)-Hy(ii+nxtot-pml_thickness-1,:,:)); % Recursive convolution
    end
    Ey(nxtot-pml_thickness+1:nxtot,:,:) = Ey(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_ey_hz_xr.*Psi_ey_hz_xr;
    Ez(nxtot-pml_thickness+1:nxtot,:,:) = Ez(nxtot-pml_thickness+1:nxtot,:,:)+CPsi_ez_hy_xr.*Psi_ez_hy_xr;
    %-------------------------------------------------------------------yleft-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_ez_hx_yl(:,jj,:) = byl(jj) * Psi_ez_hx_yl(:,jj,:)+ayl(jj)*(Hx(:,jj+1,:)-Hx(:,jj,:)); % Recursive convolution
        Psi_ex_hz_yl(:,jj,:) = byl(jj) * Psi_ex_hz_yl(:,jj,:)+ayl(jj)*(Hz(:,jj+1,:)-Hz(:,jj,:)); % Recursive convolution
    end
    Ez(:,2:pml_thickness+1,:) = Ez(:,2:pml_thickness+1,:)+CPsi_ez_hx_yl.*Psi_ez_hx_yl;
    Ex(:,2:pml_thickness+1,:) = Ex(:,2:pml_thickness+1,:)+CPsi_ex_hz_yl.*Psi_ex_hz_yl;
    %-------------------------------------------------------------------yright-----------------------------------------------------------------------------------
    for jj = 1:pml_thickness
        Psi_ez_hx_yr(:,jj,:) = byr(jj)*Psi_ez_hx_yr(:,jj,:)+ayr(jj)*(Hx(:,jj+nytot-pml_thickness,:)-Hx(:,jj+nytot-pml_thickness-1,:)); % Recursive convolution
        Psi_ex_hz_yr(:,jj,:) = byr(jj)*Psi_ex_hz_yr(:,jj,:)+ayr(jj)*(Hz(:,jj+nytot-pml_thickness,:)-Hz(:,jj+nytot-pml_thickness-1,:)); % Recursive convolution
    end
    Ez(:,nytot-pml_thickness+1:nytot,:) = Ez(:,nytot-pml_thickness+1:nytot,:)+CPsi_ez_hx_yr.*Psi_ez_hx_yr;
    Ex(:,nytot-pml_thickness+1:nytot,:) = Ex(:,nytot-pml_thickness+1:nytot,:)+CPsi_ex_hz_yr.*Psi_ex_hz_yr;
    %-------------------------------------------------------------------zleft-----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_ex_hy_zl(:,:,kk) = bzl(kk)*Psi_ex_hy_zl(:,:,kk)+azl(kk)*(Hy(:,:,kk+1)-Hy(:,:,kk)); % Recursive convolution
        Psi_ey_hx_zl(:,:,kk) = bzl(kk)*Psi_ey_hx_zl(:,:,kk)+azl(kk)*(Hx(:,:,kk+1)-Hx(:,:,kk)); % Recursive convolution
    end
    
    Ex(:,:,2:pml_thickness+1) = Ex(:,:,2:pml_thickness+1)+CPsi_ex_hy_zl.*Psi_ex_hy_zl;
    Ey(:,:,2:pml_thickness+1) = Ey(:,:,2:pml_thickness+1)+CPsi_ey_hx_zl.*Psi_ey_hx_zl;
    %-------------------------------------------------------------------zright-----------------------------------------------------------------------------------
    for kk = 1:pml_thickness
        Psi_ex_hy_zr(:,:,kk) = bzr(kk)*Psi_ex_hy_zr(:,:,kk)+azr(kk)*(Hy(:,:,kk+nztot-pml_thickness)-Hy(:,:,kk+nztot-pml_thickness-1)); % Recursive convolution
        Psi_ey_hx_zr(:,:,kk) = bzr(kk)*Psi_ey_hx_zr(:,:,kk)+azr(kk)*(Hx(:,:,kk+nztot-pml_thickness)-Hx(:,:,kk+nztot-pml_thickness-1)); % Recursive convolution
    end
    Ex(:,:,nztot-pml_thickness+1:nztot) = Ex(:,:,nztot-pml_thickness+1:nztot)+CPsi_ex_hy_zr.*Psi_ex_hy_zr;
    Ey(:,:,nztot-pml_thickness+1:nztot) = Ey(:,:,nztot-pml_thickness+1:nztot)+CPsi_ey_hx_zr.*Psi_ey_hx_zr;
    %-------------------------------------------------------------------end-----------------------------------------------------------------------------------
    
    %% Calculate Magnetic and Electric Currents
    
    %     %------------------------- xp face magnetic and electric currents , we have J_y=-H_z , J_z=H_y , M_y=E_z , M_z=-E_y--------------------------------------
    %
    %     tdmyxp(1,1,:,:) =  0.5*(Ez(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ez(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    %     tdmzxp(1,1,:,:) = -0.5*(Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %     tdjyxp(1,1,:,:) =-0.25*(Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke) ...
    %         + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %
    %     tdjzxp(1,1,:,:) =0.25*(Hy(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hy(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1) ...
    %         + Hy (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hy (ffield_ie-1,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    %
    %     %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %     %------------------------ xn face magnetic and electric currents , we have J_y=H_z , J_z=-H_y , M_y=-E_z , M_z=+E_y-------------------------------------
    %
    %     tdmyxn(1,1,:,:) = -0.5 * (Ez(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ez(ffield_is,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    %     tdmzxn(1,1,:,:) =  0.5 * (Ey(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ey(ffield_is,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %
    %     tdjyxn(1,1,:,:) = 0.25*(Hz(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hz(ffield_is,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke) ...
    %         + Hz (ffield_is-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_is-1,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %
    %     tdjzxn(1,1,:,:) =-0.25*(Hy(ffield_is,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hy(ffield_is,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1) ...
    %         + Hy (ffield_is-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hy (ffield_is-1,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    %
    %     %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %
    %     %------------------------yp face magnetic and electric currents , we have J_x=H_z , J_z=-H_x , M_x=-E_z , M_z=E_x-------------------------------------
    %
    %     tdmxyp(1,:,1,:) = -0.5*(Ez(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Ez(ffield_is+1:ffield_ie,ffield_je,ffield_ks:ffield_ke-1));
    %     tdmzyp(1,:,1,:) =  0.5*(Ex(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Ex(ffield_is:ffield_ie-1,ffield_je,ffield_ks+1:ffield_ke));
    %
    %     tdjxyp(1,:,1,:) =  0.25*(Hz(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Hz(ffield_is:ffield_ie-1,ffield_je,ffield_ks+1:ffield_ke) ...
    %         + Hz (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks+1:ffield_ke));
    %     tdjzyp(1,:,1,:) =-0.25*(Hx(ffield_is:ffield_ie-1,ffield_je,ffield_ks:ffield_ke-1)+Hx(ffield_is+1:ffield_ie,ffield_je,ffield_ks:ffield_ke-1) ...
    %         + Hx (ffield_is:ffield_ie-1,ffield_je-1,ffield_ks:ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_je-1,ffield_ks:ffield_ke-1));
    %
    %     %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %
    %     %------------------------yn face magnetic and electric currents , we have J_x=-H_z , J_z=H_x , M_x=E_z , M_z=-E_x-------------------------------------
    %     %
    %     tdmxyn(1,:,1,:) =  0.5 * (Ez(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Ez(ffield_is+1:ffield_ie,ffield_js,ffield_ks:ffield_ke-1));
    %     tdmzyn(1,:,1,:) = -0.5 * (Ex(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Ex(ffield_is:ffield_ie-1,ffield_js,ffield_ks+1:ffield_ke));
    %
    %
    %     tdjzyn(1,:,1,:) = 0.25*(Hx(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Hx(ffield_is+1:ffield_ie,ffield_js,ffield_ks:ffield_ke-1) ...
    %         + Hx (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks:ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_js-1,ffield_ks:ffield_ke-1));
    %
    %     tdjxyn(1,:,1,:) =-0.25*(Hz(ffield_is:ffield_ie-1,ffield_js,ffield_ks:ffield_ke-1)+Hz(ffield_is:ffield_ie-1,ffield_js,ffield_ks+1:ffield_ke) ...
    %         + Hz (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks:ffield_ke-1) + Hz (ffield_is:ffield_ie-1,ffield_js-1,ffield_ks+1:ffield_ke));
    %
    %     %------------------------- zp face magnetic and electric currents , we have J_x=-H_y , J_y=H_x , M_x=E_y , M_y=-E_x--------------------------------------
    %
    %     tdmxzp(1,:,:,1) =  0.5*(Ey(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Ey(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke));
    %     tdmyzp(1,:,:,1) = -0.5*(Ex(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Ex(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke));
    %
    %     tdjyzp(1,:,:,1) = 0.25*(Hx(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke) ...
    %         + Hx (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke-1) + Hx (ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ke-1));
    %
    %     tdjxzp(1,:,:,1) =-0.25*(Hy(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke)+Hy(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke) ...
    %         + Hy (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ke-1) + Hy (ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ke-1));
    %
    %     %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %     %------------------------- zn face magnetic and electric currents , we have J_x=H_y , J_y=-H_x , M_x=-E_y , M_y=E_x-------------------------------------
    %
    %     tdmxzn(1,:,:,1) = -0.5 * (Ey(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Ey(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks));
    %     tdmyzn(1,:,:,1) =  0.5 * (Ex(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Ex(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks));
    %
    %     tdjyzn(1,:,:,1) =-0.25*(Hx(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks) ...
    %         + Hx (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks-1)+Hx(ffield_is+1:ffield_ie,ffield_js:ffield_je-1,ffield_ks-1));
    %
    %     tdjxzn(1,:,:,1) = 0.25*(Hy(ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks)+Hy(ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks) ...
    %         + Hy (ffield_is:ffield_ie-1,ffield_js:ffield_je-1,ffield_ks-1) + Hy (ffield_is:ffield_ie-1,ffield_js+1:ffield_je,ffield_ks-1));
    %     %-------------------------------------------------------------------------------------------------------------------------------------------------------
    %     %-------------------------------------Transform Magnetic and Electric currents to Frequency Domain------------------------------------------------------
    %     farfield_w=freq*2*pi;
    %
    %     for ff=1:size(freq,2)
    %         kernelh = exp(-1i*farfield_w(ff)*(n-0.5)*delt)*delt; % Exponential term for DFT , here we take 0.5 time index diff. for electric currents since they
    %         % depend on magnetic fields
    %
    %         % +- x face
    %         fdjyxp(ff,:,:,:) = fdjyxp(ff,:,:,:) + kernelh * tdjyxp(1,:,:,:);
    %         fdjzxp(ff,:,:,:) = fdjzxp(ff,:,:,:) + kernelh * tdjzxp(1,:,:,:);
    %         fdjzxn(ff,:,:,:) = fdjzxn(ff,:,:,:) + kernelh * tdjzxn(1,:,:,:);
    %         fdjyxn(ff,:,:,:) = fdjyxn(ff,:,:,:) + kernelh * tdjyxn(1,:,:,:);
    %         % +- y face
    %         fdjxyp(ff,:,:,:) = fdjxyp(ff,:,:,:) + kernelh * tdjxyp(1,:,:,:);
    %         fdjzyp(ff,:,:,:) = fdjzyp(ff,:,:,:) + kernelh * tdjzyp(1,:,:,:);
    %         fdjxyn(ff,:,:,:) = fdjxyn(ff,:,:,:) + kernelh * tdjxyn(1,:,:,:);
    %         fdjzyn(ff,:,:,:) = fdjzyn(ff,:,:,:) + kernelh * tdjzyn(1,:,:,:);
    %         % +- z face
    %         fdjyzp(ff,:,:,:) = fdjyzp(ff,:,:,:) + kernelh * tdjyzp(1,:,:,:);
    %         fdjxzp(ff,:,:,:) = fdjxzp(ff,:,:,:) + kernelh * tdjxzp(1,:,:,:);
    %         fdjyzn(ff,:,:,:) = fdjyzn(ff,:,:,:) + kernelh * tdjyzn(1,:,:,:);
    %         fdjxzn(ff,:,:,:) = fdjxzn(ff,:,:,:) + kernelh * tdjxzn(1,:,:,:);
    %
    %         kernele = exp(-1i*farfield_w(ff)*n*delt)*delt;  % Exponential term for DFT , for Magnetic currents
    %
    %         % +- x face
    %         fdmxyp(ff,:,:,:) = fdmxyp(ff,:,:,:) + kernele * tdmxyp(1,:,:,:);
    %         fdmzyp(ff,:,:,:) = fdmzyp(ff,:,:,:) + kernele * tdmzyp(1,:,:,:);
    %         fdmxyn(ff,:,:,:) = fdmxyn(ff,:,:,:) + kernele * tdmxyn(1,:,:,:);
    %         fdmzyn(ff,:,:,:) = fdmzyn(ff,:,:,:) + kernele * tdmzyn(1,:,:,:);
    %         % +- y face
    %         fdmyxp(ff,:,:,:) = fdmyxp(ff,:,:,:) + kernele * tdmyxp(1,:,:,:);
    %         fdmzxp(ff,:,:,:) = fdmzxp(ff,:,:,:) + kernele * tdmzxp(1,:,:,:);
    %         fdmyxn(ff,:,:,:) = fdmyxn(ff,:,:,:) + kernele * tdmyxn(1,:,:,:);
    %         fdmzxn(ff,:,:,:) = fdmzxn(ff,:,:,:) + kernele * tdmzxn(1,:,:,:);
    %         % +- z face
    %         fdmyzp(ff,:,:,:) = fdmyzp(ff,:,:,:) + kernele * tdmyzp(1,:,:,:);
    %         fdmxzp(ff,:,:,:) = fdmxzp(ff,:,:,:) + kernele * tdmxzp(1,:,:,:);
    %         fdmxzn(ff,:,:,:) = fdmxzn(ff,:,:,:) + kernele * tdmxzn(1,:,:,:);
    %         fdmyzn(ff,:,:,:) = fdmyzn(ff,:,:,:) + kernele * tdmyzn(1,:,:,:);
    %     end
    %-------------------------------------------------------------End of Transformation------------------------------------------------------
    if(mod(n, 1) == 0)    
        subplot 121
        slice(Exincident_current,30,[],[]);  shading interp;
        caxis([-Exincamp Exincamp]);
        %colorbar; 
        axis equal;
        view([90,0]); 
        drawnow;
        F = getframe(fig);
        writeVideo(vidObj, F);
        subplot 122
        slice(Ex,30,[],[]);   shading interp;
        caxis([-0.25    0.25]);
        %colorbar; 
        axis equal;
        view([90,0]); 
        drawnow;   
        F = getframe(fig);
        writeVideo(vidObj, F);
        clc;
        disp(['Total number of computation cells (', num2str(nxtot), 'x', num2str(nytot), 'x',...
            num2str(nztot), ') : ', num2str(nxtot*nytot*nztot)]);
        disp(['FDTD progress: ',num2str((n/T)*100), ' %']);
    end
end