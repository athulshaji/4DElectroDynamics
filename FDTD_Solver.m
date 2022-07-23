function [  tdjyxp, tdjzxp, tdmyxp, tdmzxp,...
    tdjyxn, tdjzxn, tdmyxn, tdmzxn, ...
    tdjxyp, tdjzyp, tdmxyp, tdmzyp, ...
    tdjxyn, tdjzyn, tdmxyn, tdmzyn, ...
    tdjxzp, tdjyzp, tdmxzp, tdmyzp, ...
    tdjxzn, tdjyzn, tdmxzn, tdmyzn, ...
    fdjyxp, fdjzxp, fdmyxp, fdmzxp, ...
    fdjyxn, fdjzxn, fdmyxn, fdmzxn, ...
    fdjxyn, fdjzyn, fdmxyn, fdmzyn, ...
    fdjxyp, fdjzyp, fdmxyp, fdmzyp, ...
    fdjxzp, fdjyzp, fdmxzp, fdmyzp, ...
    fdjxzn, fdjyzn, fdmxzn, fdmyzn] = FDTD(     T, ...
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
    bzl, bzml, azl, azml, bzr, bzmr, azr, azmr, ...
    ffield_is, ffield_js, ffield_ks, ...
    ffield_ie, ffield_je, ffield_ke, ...
    freqs_of_interest, ...
    tdjyxp, tdjzxp, tdmyxp, tdmzxp,...
    tdjyxn, tdjzxn, tdmyxn, tdmzxn, ...
    tdjxyp, tdjzyp, tdmxyp, tdmzyp, ...
    tdjxyn, tdjzyn, tdmxyn, tdmzyn, ...
    tdjxzp, tdjyzp, tdmxzp, tdmyzp, ...
    tdjxzn, tdjyzn, tdmxzn, tdmyzn, ...
    fdjyxp, fdjzxp, fdmyxp, fdmzxp, ...
    fdjyxn, fdjzxn, fdmyxn, fdmzxn, ...
    fdjxyn, fdjzyn, fdmxyn, fdmzyn, ...
    fdjxyp, fdjzyp, fdmxyp, fdmzyp, ...
    fdjxzp, fdjyzp, fdmxzp, fdmyzp, ...
    fdjxzn, fdjyzn, fdmxzn, fdmyzn, ...
    portx, porty, portz, Zl, Z0,...
    delx, dely, delz,eps_x,eps_y,eps_z,sigma_x,sigma_y,sigma_z,freq)

% global sim_domain_materials ISSPLIT INCLUDEMETAL ISSHORT ISDIODE ISCAP freq;
global Idx_metal_region
S=1;
CIRCUIT=1;% set to 1 if a circuit is to be analysed
FREQ=freq;
dx = delx; dy = dely; dz = delz;
Nx=nxtot;Ny=nytot;Nz=nztot;
INPUT_AMP=My_dBm2Watt(0)/(delz*4);
VIDEO = 1;
if(VIDEO)
    fig = figure('Color', 'w', 'units', 'normalized', 'OuterPosition',[0 0 1 1]);
    movieName = 'FDTD_ver1_23June2022.avi';
    videoObj = VideoWriter(movieName);
    open(videoObj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For SPICE
if(1)
    time_slice=10000;
    delt_ckt=delt/time_slice;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ickt=0;index=0;V0=[];V1=[];V2=[];I_branch=0;V2pre=0;
Vnode12pre=0;
Vnode34pre=0;
Vnode40pre=0;
ILmatchpre=0;
input_source=0;
Vpre=0;Ipre=0;diode_voltage=0;
IId=[];
tol_value=1e-5;
LHSckt=[0;0;0;0;0;0;0;0;];
Vdiodepre=0;
LHScktpre=LHSckt;
[MNAckt,LHSckt,LHScktpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload]=SPICE_Init(delt_ckt);
RHSckt=LHSckt;
Rs=20;

tic
for n=1:T

    Hxincident_past = Hxincident_current;
    Hyincident_past = Hyincident_current;
    Hzincident_past = Hzincident_current;
    Exincident_past = Exincident_current;
    Eyincident_past = Eyincident_current;
    Ezincident_past = Ezincident_current;


    % Differentiated Gaussian
    %     Exincident_current = Exincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEx).*exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    %     Eyincident_current = Eyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEy).*exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    %     Ezincident_current = Ezincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt)-t0-krEz).*exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    %     Hxincident_current = Hxincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHx).*exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    %     Hyincident_current = Hyincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHy).*exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    %     Hzincident_current = Hzincamp * (-sqrt(2*exp(1))/tau)*(((n-1)*delt+delt/2)-t0-krHz).*exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);

    % Gaussian Source
    %     Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    %     Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    %     Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    %     Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    %     Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    %     Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);

    %     Sinusoidal Source

    Exincident_current = Exincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEx));
    Eyincident_current = Eyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEy));
    Ezincident_current = Ezincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEz));
    Hxincident_current = Hxincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHx));
    Hyincident_current = Hyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHy));
    Hzincident_current = Hzincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHz));


    %             if(n<=ceil(T/4))
    %                 Exincident_current = Exincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEx));
    %                 Eyincident_current = Eyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEy));
    %                 Ezincident_current = Ezincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEz));
    %                 Hxincident_current = Hxincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHx));
    %                 Hyincident_current = Hyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHy));
    %                 Hzincident_current = Hzincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHz));
    %             else
    %                 Exincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEx));
    %                 Eyincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEy));
    %                 Ezincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0-krEz));
    %                 Hxincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHx));
    %                 Hyincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHy));
    %                 Hzincident_current = 0 .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHz));
    %             end

    %     Modulated Gaussian Source
    %     Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %     Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %     Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %     Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
    %     Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
    %     Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));

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

    %         xloc=123;yloc=122;zloc=63;
    %         Eypre=Ey(xloc,yloc,zloc);

    Ex(1:nxtot,2:nytot,2:nztot) = C1Ex(1:nxtot,2:nytot,2:nztot).*Ex(1:nxtot,2:nytot,2:nztot) ...
        + C2Ex(1:nxtot,2:nytot,2:nztot).*(Hz(1:nxtot,2:nytot,2:nztot)-Hz(1:nxtot,1:nytot-1,2:nztot)) ...
        + C3Ex(1:nxtot,2:nytot,2:nztot).*(Hy(1:nxtot,2:nytot,2:nztot)-Hy(1:nxtot,2:nytot,1:nztot-1));

    Ey(2:nxtot,1:nytot,2:nztot)=C1Ey(2:nxtot,1:nytot,2:nztot).*Ey(2:nxtot,1:nytot,2:nztot) ...
        + C2Ey(2:nxtot,1:nytot,2:nztot).*(Hx(2:nxtot,1:nytot,2:nztot)-Hx(2:nxtot,1:nytot,1:nztot-1)) ...
        + C3Ey(2:nxtot,1:nytot,2:nztot).*(Hz(2:nxtot,1:nytot,2:nztot)-Hz(1:nxtot-1,1:nytot,2:nztot));

    Ez(2:nxtot,2:nytot,1:nztot)=C1Ez(2:nxtot,2:nytot,1:nztot).*Ez(2:nxtot,2:nytot,1:nztot) ...
        + C2Ez(2:nxtot,2:nytot,1:nztot).*(Hy(2:nxtot,2:nytot,1:nztot)-Hy(1:nxtot-1,2:nytot,1:nztot)) ...
        + C3Ez(2:nxtot,2:nytot,1:nztot).*(Hx(2:nxtot,2:nytot,1:nztot)-Hx(2:nxtot,1:nytot-1,1:nztot));

    YRANGE=16;XRANGE=73:79;ZRANGE=17:19;
    for yloc=YRANGE %28:32
        for xloc=XRANGE
            for zloc=ZRANGE
                Ez(xloc,yloc,zloc)=INPUT_AMP*sin(2*pi*(FREQ)*(((n-1)*delt)-t0));
            end
        end
    end

    %START OF CIRCUIT
    if(CIRCUIT) 
        CIRCUIT_INPUT_CURRENT1(n)=sum(sum(Hx(60:90,90:95,19)*delx))-sum(sum(Hx(60:90,90:95,25))*delx)+...
            sum(sum(Hz(60,90:95,19:25)*delz))-sum(sum(Hz(90,90:95,19:25)*delz));
        input_source=CIRCUIT_INPUT_CURRENT1(n);
        %         CIRCUIT_INPUT_CURRENT2(n)=sum(sum(Hx(72:80,135,19)*delx))-sum(sum(Hx(72:80,135,25))*delx)+...
        %             sum(sum(Hz(72,135,19:25)*delz))-sum(sum(Hz(80,135,19:25)*delz));
        %         input_source(2)=CIRCUIT_INPUT_CURRENT2(n);
        time=(n-1)*delt;
        [MNAckt,RHSckt,LHSckt,LHScktpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current,vd,I]=MySPICE(MNAckt,RHSckt,LHSckt,LHScktpre,Cdiode,Cload,Rload,Cmatch,Lmatch,V0,V1,V2,G0,I0,tol_value,input_source,time,delt_ckt,delt,index,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre);
        %[MNAckt,RHSckt,LHSckt,LHScktpre,V0,V1,V2,G0,I0,tol_value,diode_voltage,index]=MySPICE_NO_MathcingNetwork(MNAckt,RHSckt,LHSckt,LHScktpre,Cdiode,Ctotal,Cload,Rload,V0,V1,V2,G0,I0,tol_value,input_source,time,delt_ckt,delt,index);
        Vd(n)=vd;Idiode(n)=Vd(n)/Rs;
        %figure(1);plot(Ez(76,12:136,18));grid on;ylim([-1.1*INPUT_AMP,+1.1*INPUT_AMP]);
        for zloc=17:19
            for yloc=96
                for xloc=73:79
                    Ez(xloc,yloc,zloc)=(V0(n)/(delz*4));  
                    Ex(xloc,yloc,zloc)=0; 
                    Ey(xloc,yloc,zloc)=0; 
                end
            end
        end
        for zloc=17:19
            for yloc=116
                for xloc=73:79
                    Ez(xloc,yloc,zloc)=(V2(n)/(delz*4));
                    Ex(xloc,yloc,zloc)=0; 
                    Ey(xloc,yloc,zloc)=0; 
                end
            end
        end      
        for zloc=17:19
            for yloc=136
                for xloc=73:79
                    Ez(xloc,yloc,zloc)=(sum(sum(Hx(60:90,135,19)*delx))-sum(sum(Hx(60:90,135,25))*delx)+...
                        sum(sum(Hz(60,135,19:25)*delz))-sum(sum(Hz(90,135,19:25)*delz)))*Z0;
                    Ex(xloc,yloc,zloc)=0; 
                    Ey(xloc,yloc,zloc)=0; 
                end
            end
        end   
    end
    

    Ex = Ex + CEx_current.* Exincident_current + CEx_past.* Exincident_past;
    Ey = Ey + CEy_current.* Eyincident_current + CEy_past.* Eyincident_past;
    Ez = Ez + CEz_current.* Ezincident_current + CEz_past.* Ezincident_past;
   

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

    %------------------------- xp face magnetic and electric currents , we have J_y=-H_z , J_z=H_y , M_y=E_z , M_z=-E_y--------------------------------------

    %     tdmyxp(1,1,:,:) =  0.5*(Ez(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ez(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));
    %     tdmzxp(1,1,:,:) = -0.5*(Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Ey(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %     tdjyxp(1,1,:,:) =-0.25*(Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hz(ffield_ie,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke) ...
    %         + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hz (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks+1:ffield_ke));
    %
    %     tdjzxp(1,1,:,:) =0.25*(Hy(ffield_ie,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1)+Hy(ffield_ie,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1) ...
    %         + Hy (ffield_ie-1,ffield_js:ffield_je-1,ffield_ks:ffield_ke-1) + Hy (ffield_ie-1,ffield_js+1:ffield_je,ffield_ks:ffield_ke-1));

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
    %     farfield_w=freqs_of_interest*2*pi;
    %
    %     for ff=1:size(freqs_of_interest,2)
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
    %     Ex_monitor1(n)=Exincident_current(41,30,41);
    %     Ex_monitor2(n)=Ex(41,41,41);
    %     Ex_monitor3(n)=Ex(41,50,41);
    %     if getappdata(f,'canceling')
    %         break;
    %     end
    %     waitbar(n/T, f);
    %kk=75;
    %     arr=[27,27,56];
    %     myPlot(sim_domain_materials,Exincident_current,Ez,-Exincamp,Exincamp,n,T,arr,fig,videoObj);
    %     Nx = nxtot;
    %     Ny = nytot;
    %     Nz = nztot;
    %     el = 0;
    %     az = 90;
    %         subplot 231
    %         slice(Exincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal
    %         subplot 232
    %         slice(Eyincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal
    %         subplot 233
    %         slice(Ezincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal

    if(VIDEO)
        val = INPUT_AMP;        
        subplot 121 %234        
        slice(Ez,ceil(nytot/2)-58,ceil(nxtot/2),[ceil(nztot/2)-1,ceil(nztot/2)+1]);view(0,90); title(['iter#',num2str(n)]);%shading interp;
        caxis([-val, val]);
        colorbar;
        axis equal
        subplot 122 %235
        val = INPUT_AMP/(120*pi/sqrt(4.4));
        slice(Hz,ceil(nytot/2)-58,ceil(nxtot/2),[ceil(nztot/2)-1,ceil(nztot/2)+1]);view(0,90); title(['iter#',num2str(n)]);%shading interp;
        caxis([-val, val]);
        colorbar;axis equal;

        FRAME = getframe(fig);
        writeVideo(videoObj, FRAME)
    end        
end
toc
Vin=V1;Vout=V2;%Idiode=diode_current;Isrc=input_source;
Iin=Idiode;Iout=Idiode;
My_MATLAB_Conversion_Efficiency(Vin,Iin,Vout,Iout,delt,FREQ)
[max(Ex(:)),max(Ey(:)),max(Ez(:))]
[min(Ex(:)),min(Ey(:)),min(Ez(:))]
%save voltage_across_gap_32k.txt voltage_across_gap -ascii;
% HarmonicBalanceSolver(n,delt,freq,T,voltage_across_gap);
if(VIDEO)
    close(videoObj);
end
end

function [vrms]=My_dBm2Watt(dBm_val)
clc;close all;
format long
watt_val=10^((dBm_val/10)-3);
vrms=sqrt(50*watt_val)*sqrt(2);
end

function [MNA,LHS,LHSpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload]=SPICE_Init(delt)
% diode parameters
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
Vdiode=0;
Vdiodepre=Vdiode;
N=1.05;
Tt=1e-11;
Cj0=0.14e-12;
Vj=0.34;
M=0.4;
Fc=0.5;
Bv=2;
Ibv=1e-4;
Xti=2;
Eg=0.69;
AREA=1;%default value
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
tol_value=1e-5;
LHS=[0;0;0;0;0;0;0;0;];
LHSpre=LHS;
if Vdiode<(-10*N*Vth)
    G0=(Isat/(N*Vth))*exp(-10);
    Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
else
    G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
    Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
end
I0=Idiodepre-G0*Vdiodepre;
% load
Cload=1e-12;
Rload=500;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;
Csrc=1e-9;
Rsrc=50;

MNA=zeros(8);
MNA(2,2)=1/Rs;
MNA(2,3)=-1/Rs;
MNA(3,2)=-1/Rs;
MNA(3,3)=(1/Rs)+G0;
MNA(3,4)=-G0;
MNA(4,3)=-G0;
MNA(4,4)=G0+(1/Rload);
%MNA matrix entry for the Cmatch
MNA(5,1)=Cmatch/delt;
MNA(5,2)=-Cmatch/delt;
MNA(5,5)=-1;
MNA(2,5)=-1;
MNA(1,5)=1;
% MNA matrix entry for the Lmatch
MNA(8,2)=1;
MNA(8,8)=-Lmatch/delt;
MNA(2,8)=1;



if Vdiode<=Fc*Vj
    Cj=AREA*Cj0*(1-Vdiode/Vj).^-M;
elseif Vdiode>Fc*Vj
    Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*(Vdiode-Fc*Vj));
end

Cdiff=Tt*G0;
Cdiode=Cdiff+Cj;
% modifing the MNA matrix
MNA(3,6)=1;
MNA(4,6)=-1;
MNA(6,3)=Cdiode/delt;
MNA(6,4)=-Cdiode/delt;
MNA(6,6)=-1;

%% MAN matrix for the load section
MNA(7,4)=Cload/delt;
MNA(7,7)=-1;
MNA(4,7)=1;
% Inductor
MNA(8,2)=1;
MNA(8,8)=-Lmatch/delt;
MNA(2,8)=1;
end

function [MNA,RHS,LHS,LHSpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current,vd,I]=MySPICE(MNA,RHS,LHS,LHSpre,Cdiode,Cload,Rload,Cmatch,Lmatch,V0,V1,V2,G0,I0,tol_value,input_source,time,delt,delt_old,index,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre)
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
Vdiode=0;
Vdiodepre=Vdiode;
N=1.05;
Tt=1e-11;
Cj0=0.14e-12;
Vj=0.34;
M=0.4;
Fc=0.5;
Bv=2;
Ibv=1e-4;
Xti=2;
Eg=0.69;
AREA=1;%default value
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
Csrc=1e-9;
SRC=0;
LHS1=[];
LHS2=[];
LHS3=[];
LHS4=[];
for n=time:delt:time+(delt_old-delt)
    tol=1e9;
    SRC=input_source;
    Vnode12pre=LHS(1)-LHS(2);
    Vnode40pre=LHS(4)-0;
    ILmatchpre=LHS(8);
    Vnode34pre=LHS(3)-LHS(4);
    Vdiode=(LHS(3)-LHS(4));
    Vdiodepre=Vdiode; % Vdiodepre is the previous value of the diode voltage
    while(tol>tol_value)
        % RHS vector
        RHS=[SRC;...
            0;
            -I0;...
            I0;...
            (Cmatch/delt)*Vnode12pre;...
            (Cdiode/delt)*Vnode34pre;...
            (Cload/delt)*Vnode40pre;...
            (-Lmatch/delt)*ILmatchpre;];

        % LHS is the solution vector
        LHS=MNA\RHS;
        Vdiode=(LHS(3)-LHS(4));%Vnode34pre=LHS(3)-LHS(4);
        Vdiodepre=Vdiode;
        Vmax=N*Vth*log((ExplI/Isat)+1);
        Vbmax=N*Vth*log(ExplI/Ibv);
        Gbmax=ExplI/(N*Vth);
        % G0 and I0 are modified using applicable equations at different voltage ranges
        if Vdiodepre<(-10*N*Vth)
            if -(Bv+Vbmax)>Vdiodepre
                Idiodepre=-(ExplI+(-((Vdiodepre)+Bv)-Vbmax)*Gbmax-Ib0);
                G0=Gbmax;
            elseif ((-Bv+(308*N*Vth)>Vdiodepre) && (Vdiodepre>=-Vbmax-Bv))
                Idiodepre=(-Ibv*exp(-((Vdiodepre)+Bv)/(N*Vth)))+Ib0;
                G0=-Idiodepre/(N*Vth);
            else
                G0=(Isat/(N*Vth))*exp(-10);
                Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+(10*N*Vth)));
            end
        elseif (-10*N*Vth<=Vdiodepre && Vdiodepre<=Vmax)
            G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
            Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
        else
            Idiodepre=0;
            G0=0;
        end
        I0=Idiodepre-G0*(LHS(3)-LHS(4));
        % junction capacitance is evaluated
        if Vdiodepre<=Fc*Vj
            Cj=AREA*Cj0*(1-(Vdiodepre)/Vj).^-M;
        elseif Vdiodepre>Fc*Vj
            Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*((Vdiodepre)-Fc*Vj));
        end

        Cdiff=Tt*G0;
        Cdiode=Cdiff+Cj;
        % MNA matrix is updated
        MNA(3,3)=(1/Rs)+G0;
        MNA(3,4)=-G0;
        MNA(4,3)=-G0;
        MNA(4,4)=G0+(1/Rload);
        MNA(6,3)=Cdiode/delt;
        MNA(6,4)=-Cdiode/delt;
        
        % computing the tolerance
        %tol=norm(LHS(1:4)-LHSpre(1:4)).^2;
        tol=abs((LHS(3)-LHS(4))-(LHSpre(3)-LHSpre(4)));
        % storing the value of the LHS to LHSpre
        LHSpre=LHS;%Vdiodepre=Vdiode;
        %LHS1=[LHS1;LHS(1)];LHS2=[LHS2;LHS(2)];LHS3=[LHS3;LHS(3)];LHS4=[LHS4;LHS(4)];
    end
end
%disp(condest(MNA));
index=index+1;
V0(index)=LHS(1);
V1(index)=LHS(2);
V2(index)=LHS(4);
I(index)=LHS(6);
I_branch=SRC;
diode_current=Idiodepre;
vd=LHS(2)-LHS(3);
src_current=SRC;
end

function [MNA,LHS,LHSpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload,Ctotal]=SPICE_Init_NO_MatchingNetwork(delt)
% diode parameters
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
Vdiode=0;
Vdiodepre=Vdiode;
N=1.05;
Tt=1e-11;
Cj0=0.14e-12;
Vj=0.34;
M=0.4;
Fc=0.5;
Bv=2;
Ibv=1e-4;
Xti=2;
Eg=0.69;
AREA=1;%default value
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
tol_value=1e-5;
LHS=[0;0;0;0;0;0;0;0;];
LHSpre=LHS;
if Vdiode<(-10*N*Vth)
    G0=(Isat/(N*Vth))*exp(-10);
    Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
else
    G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
    Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
end
I0=Idiodepre-G0*Vdiodepre;
% load
Cload=1e-12;
Rload=500;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;
Csrc=1e-9;
Rsrc=50;

MNA=zeros(5);
MNA(1,1)=1/Rs;
MNA(1,2)=-1/Rs;
MNA(2,1)=-1/Rs;
MNA(2,2)=(1/Rs)+G0;
MNA(2,3)=-G0;
MNA(3,2)=-G0;
MNA(3,3)=G0+(1/Rload);

if Vdiode<=Fc*Vj
    Cj=AREA*Cj0*(1-Vdiode/Vj).^-M;
elseif Vdiode>Fc*Vj
    Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*(Vdiode-Fc*Vj));
end

Cdiff=Tt*G0;
Cdiode=Cdiff+Cj;
% modifing the MNA matrix
MNA(2,4)=1;
MNA(3,4)=-1;
MNA(4,2)=Cdiode/delt;
MNA(4,3)=-Cdiode/delt;
MNA(4,4)=-1;

%% MNA matrix for the load section
MNA(5,3)=Cload/delt;
MNA(5,5)=-1;
MNA(3,5)=1;

% Ctotal in parallel with source
Ctotal=0.1168728e-12;
MNA(6,1)=Ctotal/delt;
MNA(6,6)=-1;
MNA(1,6)=1;

Ctotal=0.1168728e-12;
MNA(7,3)=Ctotal/delt;
MNA(7,7)=-1;
MNA(3,7)=1;
end

function [MNA,RHS,LHS,LHSpre,V0,V1,V2,G0,I0,tol_value,diode_voltage,index]=MySPICE_NO_MathcingNetwork(MNA,RHS,LHS,LHSpre,Cdiode,Ctotal,Cload,Rload,V0,V1,V2,G0,I0,tol_value,input_source,time,delt,delt_old,index)
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
Vdiode=0;
Vdiodepre=Vdiode;
N=1.05;
Tt=1e-11;
Cj0=0.14e-12;
Vj=0.34;
M=0.4;
Fc=0.5;
Bv=2;
Ibv=1e-4;
Xti=2;
Eg=0.69;
AREA=1;%default value
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
Csrc=1e-9;
SRC=0;
LHS1=[];
LHS2=[];
LHS3=[];
LHS4=[];
for n=time:delt:time+(delt_old-delt)
    tol=1e9;
    SRC1=input_source(1);
    SRC2=input_source(2);
    Vnode10pre=LHS(1)-0;
    Vnode30pre=LHS(3)-0;
    Vnode23pre=LHS(2)-LHS(3);
    Vdiode=(LHS(2)-LHS(3));
    Vdiodepre=Vdiode; % Vdiodepre is the previous value of the diode voltage
    while(tol>tol_value)
        % RHS vector
        RHS=[SRC1;...
            -I0;...
            I0+SRC2;...
            (Cdiode/delt)*Vnode23pre;...
            (Cload/delt)*Vnode30pre;...
            (Ctotal/delt)*Vnode10pre;
            (Ctotal/delt)*Vnode30pre;];

        % LHS is the solution vector
        LHS=MNA\RHS;
        Vdiode=(LHS(2)-LHS(3));%Vnode34pre=LHS(3)-LHS(4);
        Vdiodepre=Vdiode;
        Vmax=N*Vth*log((ExplI/Isat)+1);
        Vbmax=N*Vth*log(ExplI/Ibv);
        Gbmax=ExplI/(N*Vth);
        % G0 and I0 are modified using applicable equations at different voltage ranges
        if Vdiodepre<(-10*N*Vth)
            if -(Bv+Vbmax)>Vdiodepre
                Idiodepre=-(ExplI+(-((Vdiodepre)+Bv)-Vbmax)*Gbmax-Ib0);
                G0=Gbmax;
            elseif ((-Bv+(308*N*Vth)>Vdiodepre) && (Vdiodepre>=-Vbmax-Bv))
                Idiodepre=(-Ibv*exp(-((Vdiodepre)+Bv)/(N*Vth)))+Ib0;
                G0=-Idiodepre/(N*Vth);
            else
                G0=(Isat/(N*Vth))*exp(-10);
                Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+(10*N*Vth)));
            end
        elseif (-10*N*Vth<=Vdiodepre && Vdiodepre<=Vmax)
            G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
            Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
        else
            Idiodepre=0;
            G0=0;
        end
        I0=Idiodepre-G0*(LHS(2)-LHS(3));
        % junction capacitance is evaluated
        if Vdiodepre<=Fc*Vj
            Cj=AREA*Cj0*(1-(Vdiodepre)/Vj).^-M;
        elseif Vdiodepre>Fc*Vj
            Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*((Vdiodepre)-Fc*Vj));
        end

        Cdiff=Tt*G0;
        Cdiode=Cdiff+Cj;
        % MNA matrix is updated
        MNA(2,2)=(1/Rs)+G0;
        MNA(2,3)=-G0;
        MNA(3,2)=-G0;
        MNA(3,3)=G0+(1/Rload);
        MNA(4,2)=Cdiode/delt;
        MNA(4,3)=-Cdiode/delt;

        % computing the tolerance
        %tol=norm(LHS(1:4)-LHSpre(1:4)).^2;
        tol=abs((LHS(2)-LHS(3))-(LHSpre(2)-LHSpre(3)))/abs((LHS(2)-LHS(3)));
        % storing the value of the LHS to LHSpre
        LHSpre=LHS;%Vdiodepre=Vdiode;
        %LHS1=[LHS1;LHS(1)];LHS2=[LHS2;LHS(2)];LHS3=[LHS3;LHS(3)];LHS4=[LHS4;LHS(4)];
    end
end
%disp(condest(MNA));
index=index+1;
V0(index)=LHS(1);
V1(index)=LHS(2);
V2(index)=LHS(3);
I(index)=SRC;
I_branch=SRC;
diode_current=Idiodepre;
diode_voltage=LHS(1)-LHS(2);
src_current=SRC;
end


function My_MATLAB_Conversion_Efficiency(Vleft,Isrc,Vright,Idiode,delt,freq)
Vleft=Vleft(1:end);
Vright=Vright(1:end);
Idiode=Idiode(1:end);
Isrc=Isrc(1:end);
NN=1e7; % appending the data with zeros for finer freq resolution
Vleft=[Vleft,zeros(1,NN)];
Idiode=[Idiode,zeros(1,NN)];
Isrc=[Isrc,zeros(1,NN)];
Vright=[Vright,zeros(1,NN)];
LL=length(Vleft);
START=1;
STOP=length(Vleft);
INPUT=Vleft(START:STOP); %(V_antenna_port/dely).^2/(2*120*pi);

%% Input
fft_Vin=fftshift(fft(INPUT))/LL;
%fft_Vin=fftshift(fft(V_txline_sample));
LL=length(fft_Vin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(1),hold on;plot(ff/1e9,log10(abs(fft_Vin)));grid on;xlim([-10,10])
fft_Iin=fftshift(fft(Isrc(START:STOP)))/LL;
LL=length(fft_Iin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(2),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pin=0.5*real(fft_Vin .* conj(fft_Iin));

%% Output
fft_Vout=fftshift(fft(Vright(START:STOP)))/LL;
LL=length(fft_Vout);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(4),hold on;plot(ff/1e9,log10(abs(fft_Vout)));grid on;xlim([-10,10])
fft_Iout=fftshift(fft(Idiode(START:STOP)))/LL;
LL=length(fft_Iout);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(5),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pout=0.5*real(fft_Vout .* conj(fft_Iout));
%figure(6),hold on;plot(ff/1e9,(abs(Pout)));grid on;xlim([-30,30])
%figure(10),hold on;plot(V_antenna_port);


% calculating the input power (RF power)
pin_val=abs(interp1(ff,(Pin),freq*0,'nearest')) +...
    2*(abs(interp1(ff,(Pin),freq*1,'nearest')) +...
    abs(interp1(ff,(Pin),freq*2,'nearest')) +...
    abs(interp1(ff,(Pin),freq*3,'nearest')) +...
    abs(interp1(ff,(Pin),freq*5,'nearest'))+...
    abs(interp1(ff,(Pin),freq*6,'nearest'))+...
    abs(interp1(ff,(Pin),freq*7,'nearest'))+...
    abs(interp1(ff,(Pin),freq*8,'nearest')));

% calculating the output power (DC power)
pout_val=abs(interp1(ff,(Pout),0,'nearest')) +...
    2*(abs(interp1(ff,(Pout),freq*1,'nearest')) +...
    abs(interp1(ff,(Pout),freq*2,'nearest')) +...
    abs(interp1(ff,(Pout),freq*3,'nearest')) +...
    abs(interp1(ff,(Pout),freq*5,'nearest'))+...
    abs(interp1(ff,(Pout),freq*6,'nearest'))+...
    abs(interp1(ff,(Pout),freq*7,'nearest'))+...
    abs(interp1(ff,(Pout),freq*8,'nearest')));

% efficiency calculation
eff=(pout_val/pin_val)*100;
disp(['RF-DC conversion efficiency: ',num2str(eff),'%']);
end


function [Z0,phase_vel,V,I,beta_1,beta_2,beta_3,r,u,q,Iprev_tele]=init_TelegrapherEquation_1D(Nz,Tsteps,R,L,G,C,Lz,RLoad,delt)
% Nz=1001;
%Tsteps=2^12;
% parameters of the tx line
%R=0;L=2500e-9;G=0;C=1e-9; % Lossless tx line
%Lz=5e-2; % 0.05 m; length of the tx line

% parameters of the source and load
%RLoad=100;
Z0=sqrt(L/C);
RS=Z0;

phase_vel=1/sqrt(L*C);
%S=0.5; % Courant factor
delz=Lz/((Nz)-1);
%delt=S*delz/phase_vel; % CFL criterion
% Initializing the voltage and current vectors
V=zeros(Nz,1);
I=zeros(Nz,1);

% Defining some important constants
beta_1=2*delt/(RS*C*delz);
beta_2=2*delt/(RLoad*C*delz);
beta_3=(2*delt)/(C*delz);
r=delt^2/(L*C*delz^2);
u=beta_3/2;
q=delt/(L*delz);
Iprev_tele=0;
end

function [V,I,Iprev,raw_data_2]=TelegrapherEquationSolver_1D(V,I,u,q,Nz,Iprev,val,delt,beta_1,beta_2,beta_3,IC)
% Voltage update
% node 1
%V(1)=(1-beta_1)*V(1)-beta_3*I(1)+beta_1*val;
V(1)=val;
% nodes from 2 to Nz-1
%     k=2:Nz-1;
%     V(k)=V(k)-u*(I(k)-I(k-1));
k=2:Nz-1;
V(k)=V(k)-u*(I(k)-I(k-1));
% last node
%raw_data3=My_LTSpiceCall_4(V(ceil(Nz/2)),delt,IC3);% for including the
% diode
%IC1=[(raw_data3.variable_mat(1,end));(raw_data3.variable_mat(2,end))];
%V1(T)=IC1(1);V2(T)=IC1(2);
%V(ceil(Nz/2))=IC3(1);
%V(ceil(Nz/2)+1)=IC3(2);
%V(ceil(Nz/2)+1)=V(ceil(Nz/2));
%k=ceil(Nz/2)+2:Nz-1;
%V(k)=V(k)-u*(I(k)-I(k-1));
%raw_data4=My_LTSpiceCall_5(V(Nz-1),delt,IC4);
%V(Nz)=raw_data4.variable_mat(1,end);
%IC2=raw_data4.variable_mat(1,end);
raw_data_2=My_LTSpiceCall_2(V(Nz-1),delt,IC);
%V(Nz)=(1-beta_2)*V(Nz)+(beta_3)*I(Nz-1);
V(Nz)=raw_data_2.variable_mat(1,end);

% Current update
k=1:Nz-1;
I(k)=I(k)-q*(V(k+1)-V(k));
Iprev=I(Nz-1);
end

function [V,I,Iprev,raw_data]=TelegrapherEquationSolver_1D_2(V,I,u,q,Nz,Iprev,val,delt,beta_1,beta_3)
% Voltage update
% node 1
V(1)=(1-beta_1)*V(1)-beta_3*I(1)+beta_1*val;
%V(1)=val;
% nodes from 2 to Nz-1
k=2:Nz-1;
V(k)=V(k)-u*(I(k)-I(k-1));
% last node
%V(Nz)=(1-beta_2)*V(Nz)+(beta_3)*I(Nz-1);
raw_data=My_LTSpiceCall_3(Iprev,delt);
V(Nz)=(raw_data.variable_mat(1,end));

% Current update
k=1:Nz-1;
I(k)=I(k)-q*(V(k+1)-V(k));
Iprev=I(Nz-1);
%time(T+1)=T*delt;
%v_storage1(T+1)=v_storage1(T)+V(Nz-2);%save 05Sept2021\\v_storage_C_1uF.txt v_storage -ascii
%     v_storage2(T+1)=V(Nz-1);
%[V(Nz-5:Nz).';I(Nz-5:end).']
%pause(1)
% plotting the result
%     subplot 211
%     plot(V,'Color','k','LineWidth',2);grid on;axis tight;
%     ylim([-2 2]);
%     title(['time step= ',num2str(T)]);
%     subplot 212
%     plot(I,'Color','k','LineWidth',2);grid on;axis tight;
%     ylim([-2 2]);
%     drawnow;
end






