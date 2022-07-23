function CN_FDTD_Solver_Rser_RC_NW(T, ...
    fmax,...
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
    Hzincident_past, delx,dely,delz,S,...
    CHx_current, CHy_current, CHz_current, ...
    CHx_past, CHy_past, CHz_past, ...
    CEx_current, CEy_current, CEz_current, ...
    CEx_past, CEy_past, CEz_past, sigma_x,sigma_y,sigma_z,eps_x,eps_y,eps_z)
% H. L. Jiang et al., "Computationally Efficient CN-PML for EM Simulations,"
% in IEEE Transactions on Microwave Theory and Techniques,
% vol. 67, no. 12, pp. 4646-4655, Dec. 2019, doi: 10.1109/TMTT.2019.2946160.

% Implementing equations 7a and 7b for CFS-PML based CN-FDTD. Background
% material conductivity and permittivity included
% Working!!!

%% Initializing MATLAB

Z0=50;
INPUT_AMP=My_dBm2Watt(0,Z0)/(delz*2);
FREQ = fmax;
VIDEO = 0; % numerals 0 or 1
PML_ENABLED = 1;
file_name = '26June2022_ver1-S=2_logscale.avi';
%% Basic parameter

c0 = 299792458; %m/s
Nx = nxtot;
Ny = nytot;
Nz = nztot;

Nxyz = Nx*Ny*Nz;
Nxy = Nx*Ny;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YRANGE=76;XRANGE=65:67;ZRANGE=36:38;
% PEC_WALL=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
% eps_x(XRANGE,YRANGE,ZRANGE)=1e100;
% sigma_x(XRANGE,YRANGE,ZRANGE)=1e100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4*pi*1e-7; % A/m
epsilon0 = 8.854e-12; % F/m
epsilon = epsilon0.*eps_x(1:Nx,1:Ny,1:Nz);
sigma_x = sigma_x(1:Nx,1:Ny,1:Nz);

NxOriginal = Nx - 2*pml_thickness;
NyOriginal = Ny - 2*pml_thickness;
NzOriginal = Nz - 2*pml_thickness;

NxPML = pml_thickness;
NyPML = NxPML;
NzPML = NxPML;

dxyz = sqrt((1/delx)^2 + (1/dely)^2 + (1/delz)^2);
%% dt computation

dtCFL = 1/(c0*dxyz);
dt = S * dtCFL;

%% Basic parameters continued
%
Nxyz = Nx*Ny*Nz;
Nxy = Nx*Ny;

ExInc = zeros(Ny,Nx,Nz);
EyInc = ExInc;
EzInc = ExInc;
HxInc = ExInc;
HyInc = ExInc;
HzInc = ExInc;

%%
krEx = krEx(1:Nx,1:Ny,1:Nz);
krEy = krEy(1:Nx,1:Ny,1:Nz);
krEz = krEz(1:Nx,1:Ny,1:Nz);
krHx = krHx(1:Nx,1:Ny,1:Nz);
krHy = krHy(1:Nx,1:Ny,1:Nz);
krHz = krHz(1:Nx,1:Ny,1:Nz);

Exincident_current = Exincident_current(1:Nx,1:Ny,1:Nz);
Eyincident_current = Eyincident_current(1:Nx,1:Ny,1:Nz);
Ezincident_current = Ezincident_current(1:Nx,1:Ny,1:Nz);
Hxincident_current = Hxincident_current(1:Nx,1:Ny,1:Nz);
Hyincident_current = Hyincident_current(1:Nx,1:Ny,1:Nz);
Hzincident_current = Hzincident_current(1:Nx,1:Ny,1:Nz);

Exincident_past = Exincident_past(1:Nx,1:Ny,1:Nz);
Eyincident_past = Eyincident_past(1:Nx,1:Ny,1:Nz);
Ezincident_past = Ezincident_past(1:Nx,1:Ny,1:Nz);
Hxincident_past = Hxincident_past(1:Nx,1:Ny,1:Nz);
Hyincident_past = Hyincident_past(1:Nx,1:Ny,1:Nz);
Hzincident_past = Hzincident_past(1:Nx,1:Ny,1:Nz);

CHx_current = CHx_current(1:Nx,1:Ny,1:Nz);
CHy_current = CHy_current(1:Nx,1:Ny,1:Nz);
CHz_current = CHz_current(1:Nx,1:Ny,1:Nz);
CEx_current = CEx_current(1:Nx,1:Ny,1:Nz);
CEy_current = CEy_current(1:Nx,1:Ny,1:Nz);
CEz_current = CEz_current(1:Nx,1:Ny,1:Nz);

CHx_past = CHx_past(1:Nx,1:Ny,1:Nz);
CHy_past = CHy_past(1:Nx,1:Ny,1:Nz);
CHz_past = CHz_past(1:Nx,1:Ny,1:Nz);
CEx_past = CEx_past(1:Nx,1:Ny,1:Nz);
CEy_past = CEy_past(1:Nx,1:Ny,1:Nz);
CEz_past = CEz_past(1:Nx,1:Ny,1:Nz);

sigma_x_vec = [];
sigma_y_vec = [];
sigma_z_vec = [];
epsilon_vec = [];
for k = 1:Nz
    sigma_x_vec = [sigma_x_vec;reshape(sigma_x(1:Nx,1:Ny,k).',Nxy,1)];
    sigma_y_vec = [sigma_y_vec;reshape(sigma_y(1:Nx,1:Ny,k).',Nxy,1)];
    sigma_z_vec = [sigma_z_vec;reshape(sigma_z(1:Nx,1:Ny,k).',Nxy,1)];
    epsilon_vec = [epsilon_vec;reshape(epsilon(1:Nx,1:Ny,k).',Nxy,1)];
end

%% Defining the PML parameters
if(PML_ENABLED)
    m = 3; R = exp(-32);
    sigmaXOPT = -((m+1)*log(R)/(2*120*pi*NxPML*delx));
    sigmaXMAX = sigmaXOPT;
    sigmaYOPT = -((m+1)*log(R)/(2*120*pi*NyPML*dely));
    sigmaYMAX = sigmaYOPT;
    sigmaZOPT = -((m+1)*log(R)/(2*120*pi*NzPML*delz));
    sigmaZMAX = sigmaZOPT;
    alphaXMax = 0.2; alphaYMax = alphaXMax; alphaZMax = alphaXMax;
    kappaXMax = 1; kappaYMax = kappaXMax; kappaZMax = kappaXMax;
    %sigmaXMax = 100; sigmaYMax = sigmaXMax; sigmaZMax = sigmaXMax;
    kappaX = ones(Ny,Nx,Nz); kappaY = kappaX; kappaZ = kappaX;

    for k = 1:Nz
        for j =1:Ny
            for i = 1:NxPML
                %alphaX(j,i,k) = alphaXMax * (i/NxPML)^m;
                alphaX(j,i,k) = alphaXMax * ((NxPML-i+1)/NxPML);
                kappaX(j,i,k) = 1 + ((kappaXMax-1) * ((NxPML-i+1)/NxPML)^m);
                sigmaX(j,i,k) = sigmaXMAX * ((NxPML-i+1)/NxPML)^m;
                PMLX(j,i,k) = 1;
            end
        end
    end
    alphaX(:,Nx-NxPML+1:Nx,1:Nz) = alphaX(:,NxPML:-1:1,1:Nz);
    kappaX(:,Nx-NxPML+1:Nx,1:Nz) = kappaX(:,NxPML:-1:1,1:Nz);
    sigmaX(:,Nx-NxPML+1:Nx,1:Nz) = sigmaX(:,NxPML:-1:1,1:Nz);
    PMLX(:,Nx-NxPML+1:Nx,1:Nz) = PMLX(:,NxPML:-1:1,1:Nz);

    for k = 1:Nz
        for i = 1:Nx
            for j = 1:NyPML
                %alphaY(j,i,k) = alphaYMax * (j/NyPML)^m;
                alphaY(j,i,k) = alphaYMax * ((NyPML - j + 1)/NyPML);
                kappaY(j,i,k) = 1+((kappaYMax-1) * ((NyPML - j + 1)/NyPML)^m);
                sigmaY(j,i,k) = sigmaYMAX * ((NyPML - j + 1)/NyPML)^m;
                PMLY(j,i,k) = 1;
            end
        end
    end
    alphaY(Ny-NyPML+1:Ny,:,:) = alphaY(NyPML:-1:1,:,:);
    kappaY(Ny-NyPML+1:Ny,:,:) = kappaY(NyPML:-1:1,:,:);
    sigmaY(Ny-NyPML+1:Ny,:,:) = sigmaY(NyPML:-1:1,:,:);
    PMLY(Ny-NyPML+1:Ny,:,:) = PMLY(NyPML:-1:1,:,:);

    for i = 1:Nx
        for j = 1:Ny
            for k = 1:NzPML
                %alphaZ(j,i,k) = alphaZMax * (k/NzPML)^m;
                alphaZ(j,i,k) = alphaZMax * ((NzPML - k + 1)/NzPML);
                kappaZ(j,i,k) = 1 + ((kappaZMax-1) * ((NzPML - k + 1)/NzPML)^m);
                sigmaZ(j,i,k) = sigmaZMAX * ((NzPML - k + 1)/NzPML)^m;
                PMLZ(j,i,k) = 1;
            end
        end
    end
    alphaZ(:,:,Nz-NzPML+1:Nz) = alphaZ(:,:,NzPML:-1:1);
    sigmaZ(:,:,Nz-NzPML+1:Nz) = sigmaZ(:,:,NzPML:-1:1);
    kappaZ(:,:,Nz-NzPML+1:Nz) = kappaZ(:,:,NzPML:-1:1);
    PMLZ(:,:,Nz-NzPML+1:Nz) = PMLZ(:,:,NzPML:-1:1);

    vX = (kappaX.*(2*epsilon0 - alphaX*dt) - sigmaX*dt)./(kappaX.*(2*epsilon0 + alphaX*dt) + sigmaX*dt);
    vY = (kappaY.*(2*epsilon0 - alphaY*dt) - sigmaY*dt)./(kappaY.*(2*epsilon0 + alphaY*dt) + sigmaY*dt);
    vZ = (kappaZ.*(2*epsilon0 - alphaZ*dt) - sigmaZ*dt)./(kappaZ.*(2*epsilon0 + alphaZ*dt) + sigmaZ*dt);

    uX = ((2*epsilon0) - (alphaX*dt))./((2*epsilon0) + (alphaX*dt));
    uY = ((2*epsilon0) - (alphaY*dt))./((2*epsilon0) + (alphaY*dt));
    uZ = ((2*epsilon0) - (alphaZ*dt))./((2*epsilon0) + (alphaZ*dt));

    rX = ((vX - uX));
    rY = ((vY - uY));
    rZ = ((vZ - uZ));

    wX = ((kappaX.^2.*((2*epsilon0) + (alphaX*dt))) + (kappaX.*sigmaX.*dt)) ./ (kappaX.*((2*epsilon0) + (alphaX*dt)));
    wY = ((kappaY.^2.*((2*epsilon0) + (alphaY*dt))) + (kappaY.*sigmaY.*dt)) ./ (kappaY.*((2*epsilon0) + (alphaY*dt)));
    wZ = ((kappaZ.^2.*((2*epsilon0) + (alphaZ*dt))) + (kappaZ.*sigmaZ.*dt)) ./ (kappaZ.*((2*epsilon0) + (alphaZ*dt)));

    RX = []; RY = []; RZ = [];
    VX = []; VY = []; VZ = [];
    UX = []; UY = []; UZ = [];
    WX = []; WY = []; WZ = [];
    XPML = []; YPML = []; ZPML = [];

    ExIncVec = [];
    EyIncVec = [];
    EzIncVec = [];
    HxIncVec = [];
    HyIncVec = [];
    HzIncVec = [];

    for k = 1:Nz
        RX = [RX;reshape(rX(:,:,k).',Nxy,1)];
        RY = [RY;reshape(rY(:,:,k).',Nxy,1)];
        RZ = [RZ;reshape(rZ(:,:,k).',Nxy,1)];

        VX = [VX;reshape(vX(:,:,k).',Nxy,1)];
        VY = [VY;reshape(vY(:,:,k).',Nxy,1)];
        VZ = [VZ;reshape(vZ(:,:,k).',Nxy,1)];

        UX = [UX;reshape(uX(:,:,k).',Nxy,1)];
        UY = [UY;reshape(uY(:,:,k).',Nxy,1)];
        UZ = [UZ;reshape(uZ(:,:,k).',Nxy,1)];

        WX = [WX;reshape(wX(:,:,k).',Nxy,1)];
        WY = [WY;reshape(wY(:,:,k).',Nxy,1)];
        WZ = [WZ;reshape(wZ(:,:,k).',Nxy,1)];

        XPML = [XPML;reshape(PMLX(:,:,k).',Nxy,1)];
        YPML = [YPML;reshape(PMLY(:,:,k).',Nxy,1)];
        ZPML = [ZPML;reshape(PMLZ(:,:,k).',Nxy,1)];
    end
    clear rX rY rZ uX uY uZ vX vY vZ
else
    WX = ones(Nxyz,1); WY = WX; WZ = WX;
end

[Dxe,Dxxe,Dxh,Dye,Dyye,Dyh,Dze,Dzze,Dzh] = DerivativeOperator(delx,dely,delz,Nx,Ny,Nz,Nxy,Nxyz);
a0 = (2.*epsilon_vec - sigma_x_vec.*dt)./(2*epsilon_vec + sigma_x_vec.*dt); % 1
a1 = dt./(2.*epsilon_vec+(sigma_x_vec*delt));
a2 = ones(Nxyz,1).*dt./(2*mu0);

tic
Dxew = Dxe.*(a2./WX);
% Dxxew = (Dxxe(i,:)./(WX.'));
Dxhw = Dxh.*(a1./WX);

Dyew = Dye.*(a2./WY);
% Dyyew = (Dyye(i,:)./(WY.'));
Dyhw = Dyh.*(a1./WY);

Dzew = Dze.*(a2./WZ);
% Dzzew = (Dzze(i,:)./(WZ.'));
Dzhw = Dzh.*(a1./WZ);
toc

I = speye(Nxyz,Nxyz);

%% Solving the set of equations

ExOld = sparse(Nxyz,1);
EyOld = sparse(Nxyz,1);
EzOld = sparse(Nxyz,1);
HxOld = sparse(Nxyz,1);
HyOld = sparse(Nxyz,1);
HzOld = sparse(Nxyz,1);

ExNew = sparse(Nxyz,1);
EyNew = sparse(Nxyz,1);
EzNew = sparse(Nxyz,1);
HxNew = sparse(Nxyz,1);
HyNew = sparse(Nxyz,1);
HzNew = sparse(Nxyz,1);

if(PML_ENABLED)
    psiEXY = sparse(Nxyz,1);
    psiEXZ = sparse(Nxyz,1);
    psiEYZ = sparse(Nxyz,1);
    psiEYX = sparse(Nxyz,1);
    psiEZY = sparse(Nxyz,1);
    psiEZX = sparse(Nxyz,1);

    psiHXZ = sparse(Nxyz,1);
    psiHXY = sparse(Nxyz,1);
    psiHYX = sparse(Nxyz,1);
    psiHYZ = sparse(Nxyz,1);
    psiHZY = sparse(Nxyz,1);
    psiHZX = sparse(Nxyz,1);
end

PHI_OLD = [ExOld; EyOld; EzOld; HxOld; HyOld; HzOld];
PHI_NEW = [ExNew; EyNew; EzNew; HxNew; HyNew; HzNew];
%A = speye(6*Nxyz,6*Nxyz);
A = spdiags([a0; a0; a0; ones(Nxyz,1); ones(Nxyz,1); ones(Nxyz,1)], 0,6*Nxyz, 6*Nxyz);
Z = sparse(Nxyz,Nxyz);

%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE WALL
YRANGE=36;XRANGE=65:67;ZRANGE=36:38;
SOURCE_WALL=[];
for yloc=YRANGE %28:32
    for xloc=XRANGE
        for zloc=ZRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            SOURCE_WALL=[SOURCE_WALL;xyzloc];
            Dyhw(xyzloc,:)=0;
            Dxhw(xyzloc,:)=0;
            Dzhw(xyzloc,:)=0;
            A(2*Nxyz+xyzloc,2*Nxyz+xyzloc)=1;
            A(0*Nxyz+xyzloc,0*Nxyz+xyzloc)=0;
            A(1*Nxyz+xyzloc,1*Nxyz+xyzloc)=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% END VOLTAGE WALL
END_VOLTAGE_WALL=[];
YRANGE=96;XRANGE=65:67;ZRANGE=37;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            END_VOLTAGE_WALL=[END_VOLTAGE_WALL;xyzloc];
            Dyhw(xyzloc,:)=0;
            Dxhw(xyzloc,:)=0;
            Dzhw(xyzloc,:)=0;
            A(2*Nxyz+xyzloc,2*Nxyz+xyzloc)=1;
            A(0*Nxyz+xyzloc,0*Nxyz+xyzloc)=0;
            A(1*Nxyz+xyzloc,1*Nxyz+xyzloc)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% VOLTAGE WALL1
VOLTAGE_WALL1=[];
YRANGE=76;XRANGE=65:67;ZRANGE=36:38;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            VOLTAGE_WALL1=[VOLTAGE_WALL1;xyzloc];
            Dyhw(xyzloc,:)=0;
            Dxhw(xyzloc,:)=0;
            Dzhw(xyzloc,:)=0;
            A(2*Nxyz+xyzloc,2*Nxyz+xyzloc)=1;
            A(0*Nxyz+xyzloc,0*Nxyz+xyzloc)=0;
            A(1*Nxyz+xyzloc,1*Nxyz+xyzloc)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%
% VOLTAGE WALL2
VOLTAGE_WALL2=[];
YRANGE=86;XRANGE=65:67;ZRANGE=37;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            VOLTAGE_WALL2=[VOLTAGE_WALL2;xyzloc];
            Dyhw(xyzloc,:)=0;
            Dxhw(xyzloc,:)=0;
            Dzhw(xyzloc,:)=0;
            A(2*Nxyz+xyzloc,2*Nxyz+xyzloc)=1;
            A(0*Nxyz+xyzloc,0*Nxyz+xyzloc)=0;
            A(1*Nxyz+xyzloc,1*Nxyz+xyzloc)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 2
CURRENT_LINE21=[];
YRANGE=95;XRANGE=64:68;ZRANGE=39;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            CURRENT_LINE21=[CURRENT_LINE21;xyzloc];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 2
CURRENT_LINE22=[];
YRANGE=95;XRANGE=64:68;ZRANGE=37;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            CURRENT_LINE22=[CURRENT_LINE22;xyzloc];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 2
CURRENT_LINE23=[];
YRANGE=95;XRANGE=64;ZRANGE=37:39;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            CURRENT_LINE23=[CURRENT_LINE23;xyzloc];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 2
CURRENT_LINE24=[];
YRANGE=95;XRANGE=68;ZRANGE=37:39;
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            CURRENT_LINE24=[CURRENT_LINE24;xyzloc];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 1
YRANGE=75;XRANGE=64:68;ZRANGE=38;
CURRENT_LINE_11=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
% YRANGE=76;XRANGE=40:66;ZRANGE=39;
% CURRENT_LINE_111=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 1
YRANGE=75;XRANGE=64:68;ZRANGE=36;
CURRENT_LINE_12=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
% YRANGE=76;XRANGE=40:66;ZRANGE=37;
% CURRENT_LINE_122=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 1
YRANGE=75;XRANGE=64;ZRANGE=36:38;
CURRENT_LINE_13=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
% YRANGE=76;XRANGE=64;ZRANGE=37:39;
% CURRENT_LINE_133=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT 1
YRANGE=75;XRANGE=68;ZRANGE=36:38;
CURRENT_LINE_14=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
% YRANGE=76;XRANGE=68;ZRANGE=37:39;
% CURRENT_LINE_144=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny);
%%%%%%%%%%%%%%%%%%%%%%%

D1 = [Z Z Z Z -Dzhw Z;...
    Z Z Z Z Z -Dxhw;...
    Z Z Z -Dyhw Z Z;...
    Z Z -Dyew Z Z Z;...
    -Dzew Z Z Z Z Z;...
    Z -Dxew Z Z Z Z;];
D1 = sparse(D1);
D2 = [  Z Z Z Z Z Dyhw;...
    Z Z Z Dzhw Z Z;...
    Z Z Z Z Dxhw Z;...
    Z Dzew Z Z Z Z;...
    Z Z Dxew Z Z Z;...
    Dyew Z Z Z Z Z;];
D2 = sparse(D2);

if(PML_ENABLED)
    PSI = [ RY.*psiEXY - RZ.*psiEXZ;...
        RZ.*psiEYZ - RX.*psiEYX;...
        RX.*psiEZX - RY.*psiEZY;...
        RZ.*psiHXZ - RY.*psiHXY;...
        RX.*psiHYX - RZ.*psiHYZ;...
        RY.*psiHZY - RX.*psiHZX;];
end
EYE = speye(6*Nxyz,6*Nxyz);

tic;
[L_STAR, U_STAR, P_STAR] = lu(EYE - D1);
[L,U,P] = lu(EYE - D2);
toc;

% tic;
% [L_STAR, U_STAR, P_STAR] = ilu(EYE - D1);
% [L,U,P] = ilu(EYE - D2);
% toc;

Z = sparse(Nxyz,1);
%% Output

if(VIDEO)
    vidObj = VideoWriter(file_name);
    open(vidObj);
end

% INDEX = 0;
% PEC_LOC =


if(VIDEO)
    fig = figure('color','w','units','normalized','outerposition',[0 0 1 1]);
end

FIELD_VALUE=[];
fieldEx=zeros(Nxyz,1);fieldEy=zeros(Nxyz,1);fieldEz=zeros(Nxyz,1);
fieldHx=zeros(Nxyz,1);fieldHy=zeros(Nxyz,1);fieldHz=zeros(Nxyz,1);
% INPUT_FIELD=[fieldEx;fieldEy;fieldEz;fieldHx;fieldHy;fieldHz;];
INPUT_FIELD=sparse(6*Nxyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For SPICE
if(1)
    time_slice=1000;
    delt_ckt=delt/time_slice;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IC=[0;0;0;0;];
Ickt=0;index=0;V0=[];V1=[];V2=[];I_branch=0;V2pre=0;
Vnode12pre=0;
Vnode34pre=0;
Vnode40pre=0;
ILmatchpre=0;
input_source=0;
Vpre=0;Ipre=0;diode_voltage=0;
IId=[];
tol_value=1e-5;
Vdiodepre=0;
%LHScktpre=LHSckt;
%[MNAckt,LHSckt,LHScktpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload,Csrc]=SPICE_Init(delt_ckt,delx,dely,delz);
[MNAckt,RHSckt,LHSckt,Rload,Rser,Cload,Csrc]=SPICE_INIT_Rser_RC_NW(delt_ckt,delx,dely,delz);
RHSckt=LHSckt;
Rs=20;

tic
for n = 1:T
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

    %     Gaussian Source
    %     Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    %     Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    %     Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    %     Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    %     Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    %     Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);

    % Modulated Gaussian Source
    %         Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %         Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %         Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
    %         Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
    %         Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
    %         Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));


    % Sinusoidal Source
    %     Exincident_current = Exincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEx));
    %     Eyincident_current = Eyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEy));
    %     Ezincident_current = Ezincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEz));
    %     Hxincident_current = Hxincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHx));
    %     Hyincident_current = Hyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHy));
    %     Hzincident_current = Hzincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHz));

    ExInc = CEx_current .* Exincident_current + CEx_past .* Exincident_past;
    EyInc = CEy_current .* Eyincident_current + CEy_past .* Eyincident_past;
    EzInc = CEz_current .* Ezincident_current + CEz_past .* Ezincident_past;
    HxInc = CHx_current .* Hxincident_current + CHx_past .* Hxincident_past;
    HyInc = CHy_current .* Hyincident_current + CHy_past .* Hyincident_past;
    HzInc = CHz_current .* Hzincident_current + CHz_past .* Hzincident_past;
    ExIncVec = [];
    EyIncVec = [];
    EzIncVec = [];
    HxIncVec = [];
    HyIncVec = [];
    HzIncVec = [];

    for k = 1:Nz
        ExIncVec = [ExIncVec;reshape(ExInc(:,:,k).',Nxy,1)];
        EyIncVec = [EyIncVec;reshape(EyInc(:,:,k).',Nxy,1)];
        EzIncVec = [EzIncVec;reshape(EzInc(:,:,k).',Nxy,1)];

        HxIncVec = [HxIncVec;reshape(HxInc(:,:,k).',Nxy,1)];
        HyIncVec = [HyIncVec;reshape(HyInc(:,:,k).',Nxy,1)];
        HzIncVec = [HzIncVec;reshape(HzInc(:,:,k).',Nxy,1)];
    end

    INC_FIELD = [ExIncVec;EyIncVec;EzIncVec;HxIncVec;HyIncVec;HzIncVec;];

    IO_FIELD=[fieldEx;fieldEy;fieldEz;fieldHx;fieldHy;fieldHz;];   

    RHS_STAR = (A+D1+(2*D2))*PHI_OLD + (PSI+INC_FIELD);
    PHI_STAR = U_STAR \ (L_STAR \ (P_STAR * RHS_STAR));
    %[PHI_STAR,flag1]=gmres(EYE-D1,RHS_STAR,10,[],size(RHS_STAR,1),L_STAR,U_STAR);    
    RHS = PHI_STAR - (D2*(PHI_OLD));
    PHI_NEW = U \ (L \ (P * (RHS)));
    %[PHI_NEW,flag2]=gmres(EYE-D2,RHS,10,[],size(RHS,1),L,U);

    

    ExNew = PHI_NEW(0*Nxyz+1:1*Nxyz,1);
    EyNew = PHI_NEW(1*Nxyz+1:2*Nxyz,1);
    EzNew = PHI_NEW(2*Nxyz+1:3*Nxyz,1);
    HxNew = PHI_NEW(3*Nxyz+1:4*Nxyz,1);
    HyNew = PHI_NEW(4*Nxyz+1:5*Nxyz,1);
    HzNew = PHI_NEW(5*Nxyz+1:6*Nxyz,1);

    
    EzNew(SOURCE_WALL)=INPUT_AMP*sin(2*pi*(FREQ)*(((n-1)*delt)));  
    %ExNew(SOURCE_WALL)=0; EyNew(SOURCE_WALL)=0;

    time=(n-1)*delt;
    CIRCUIT_CURRENT_INPUT(n)=((sum(HxNew(CURRENT_LINE_11))*delx)-(sum(HxNew(CURRENT_LINE_12))*delx)+(sum(HxNew(CURRENT_LINE_13))*delz)-(sum(HxNew(CURRENT_LINE_14))*delz));
    %     CIRCUIT_CURRENT_INPUT(n)=((sum(HxNew(CURRENT_LINE_11))*delx)-(sum(HxNew(CURRENT_LINE_111))*delx))-...
    %         ((sum(HxNew(CURRENT_LINE_12))*delx)-(sum(HxNew(CURRENT_LINE_122))*delx))-...
    %         (sum(HzNew(CURRENT_LINE_13))*delz)-...
    %         (sum(HzNew(CURRENT_LINE_14))*delz);
    input_source=CIRCUIT_CURRENT_INPUT(n);
    %     [MNAckt,RHSckt,LHSckt,LHScktpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current,vd,I]=MySPICE(MNAckt,RHSckt,LHSckt,LHScktpre,Cdiode,Cload,Rload,Cmatch,Lmatch,Csrc,V0,V1,V2,G0,I0,tol_value,input_source,time,delt_ckt,delt,index,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre);
    %     Vd(n)=vd;Idiode(n)=vd/Rs;

    %[MNAckt,RHSckt,LHSckt]=MySPICE2(MNAckt,RHSckt,LHSckt,V1,V2,input_source,delt_ckt,delt,Rload,Rser,Cload,Csrc,time);
    V0(n)=Z0*input_source;  %LHSckt(1);
    V2(n)=Z0*input_source;  %LHSckt(2);
    EzNew(VOLTAGE_WALL1)=(V0(n)/(delz*2)); %ExNew(VOLTAGE_WALL1)=0; EyNew(VOLTAGE_WALL1)=0;
    EzNew(VOLTAGE_WALL2)=(V2(n)/(delz*2)); %ExNew(VOLTAGE_WALL2)=0; EyNew(VOLTAGE_WALL2)=0;

    %EzNew(PEC_WALL1)=0;  ExNew(PEC_WALL1)=0; EyNew(PEC_WALL1)=0;   

    EzNew(END_VOLTAGE_WALL)=Z0*((sum(HxNew(CURRENT_LINE21))*delx)-(sum(HxNew(CURRENT_LINE22))*delx)+(sum(HzNew(CURRENT_LINE23))*delz)-(sum(HzNew(CURRENT_LINE24))*delz));
    %ExNew(END_VOLTAGE_WALL)=0; EyNew(END_VOLTAGE_WALL)=0;    

    ExOld = PHI_OLD(0*Nxyz+1:1*Nxyz,1);
    EyOld = PHI_OLD(1*Nxyz+1:2*Nxyz,1);
    EzOld = PHI_OLD(2*Nxyz+1:3*Nxyz,1);
    HxOld = PHI_OLD(3*Nxyz+1:4*Nxyz,1);
    HyOld = PHI_OLD(4*Nxyz+1:5*Nxyz,1);
    HzOld = PHI_OLD(5*Nxyz+1:6*Nxyz,1);

    psiEXY = VY.*psiEXY + Dyhw*((HzNew + HzOld));
    psiEXZ = VZ.*psiEXZ + Dzhw*((HyNew + HyOld));
    psiEYZ = VZ.*psiEYZ + Dzhw*((HxNew + HxOld));
    psiEYX = VX.*psiEYX + Dxhw*((HzNew + HzOld));
    psiEZX = VX.*psiEZX + Dxhw*((HyNew + HyOld));
    psiEZY = VY.*psiEZY + Dyhw*((HxNew + HxOld));

    psiHXY = VY.*psiHXY + Dyew*((EzNew + EzOld));
    psiHXZ = VZ.*psiHXZ + Dzew*((EyNew + EyOld));
    psiHYZ = VZ.*psiHYZ + Dzew*((ExNew + ExOld));
    psiHYX = VX.*psiHYX + Dxew*((EzNew + EzOld));
    psiHZX = VX.*psiHZX + Dxew*((EyNew + EyOld));
    psiHZY = VY.*psiHZY + Dyew*((ExNew + ExOld));

    PHI_OLD = [ExNew;EyNew;EzNew;HxNew;HyNew;HzNew;];   

    PSI = [ RY.*psiEXY - RZ.*psiEXZ ;...
        RZ.*psiEYZ - RX.*psiEYX ;...
        RX.*psiEZX - RY.*psiEZY ;...
        RZ.*psiHXZ - RY.*psiHXY ;...
        RX.*psiHYX - RZ.*psiHYZ ;...
        RY.*psiHZY - RX.*psiHZX ;];


    for i = 1:Nz
        FIELDEX(:,:,i) = full(reshape(ExNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');
        FIELDHX(:,:,i) = full(reshape(HxNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');

        FIELDEY(:,:,i) = full(reshape(EyNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');
        FIELDHY(:,:,i) = full(reshape(HyNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');

        FIELDEZ(:,:,i) = full(reshape(EzNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');
        FIELDHZ(:,:,i) = full(reshape(HzNew((i-1)*Nxy+1:i*Nxy,1),Ny,Nx).');
    end
%     field_probe=FIELDEZ(91,:,33);
%     figure(2),grid on;plot(field_probe);ylim([-2.1*INPUT_AMP,2.1*INPUT_AMP]);title(['iter# ',num2str(n)]);


    if(VIDEO && n>3)
        val = log10((INPUT_AMP));
        %subplot 121 %234
        slice(log10(abs(FIELDEZ)),ceil(nytot/2),ceil(nxtot/2),[ceil(nztot/2)-1,ceil(nztot/2)+2]);view(0,90); title(['iter#',num2str(n)]);%shading interp;
        caxis([-val, val]);
        colorbar;
        axis equal
        %         subplot 122 %235
        %         val = INPUT_AMP/(120*pi/sqrt(4.4));
        %         slice(FIELDHZ,ceil(nytot/2),ceil(nxtot/2),[ceil(nztot/2)-1,ceil(nztot/2)+1]);view(0,90); title(['iter#',num2str(n)]);%shading interp;
        %         caxis([-val, val]);
        %         colorbar;axis equal;

        FRAME = getframe(fig);
        writeVideo(vidObj, FRAME)
    end
    clc;disp(n);
end
toc
if(VIDEO)
    close(vidObj);
end
Vin=V1;Vout=V2;
Iin=Idiode;Iout=Idiode;
My_MATLAB_Conversion_Efficiency(Vin,Iin,Vout,Iout,delt,FREQ)
Vin=V0;Vout=V2;
Iin=CIRCUIT_CURRENT_INPUT;Iout=Idiode;
My_MATLAB_Conversion_Efficiency(Vin,Iin,Vout,Iout,delt,FREQ)
figure,hold on;grid on;plot(0:delt:(n-1)*delt,V1,'r');plot(0:delt:(n-1)*delt,V2,'k');drawnow;
end

function PATH=ThreeDimNodeLoc2VectorNodesLoc(XRANGE,YRANGE,ZRANGE,Nxy,Ny)
PATH=[];
for zloc=ZRANGE
    for yloc=YRANGE
        for xloc=XRANGE
            xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
            PATH=[PATH;xyzloc];
        end
    end
end
end



function [Vpeak]=My_dBm2Watt(dBm_val,Z0)
clc;close all;
format long
watt_val=10^((dBm_val/10)-3);
vrms=sqrt(Z0*watt_val);
Vpeak=vrms*sqrt(2);
end

function [MNA,LHS,LHSpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload,Csrc]=SPICE_Init(delt,delx,dely,delz)
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
LHS=[0;0;0;0;0;0;0;0;0;];
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
cap_area=dely*delx;
cap_height=delz*2;
Csrc=8.854e-12*4.4*(cap_area)/(cap_height);
% Csrc=3*8.854e-12*4.4*0.002;
Rsrc=50;

MNA=zeros(9);
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
%% MAN matrix for the load section
MNA(9,1)=Csrc/delt;
MNA(9,9)=-1;
MNA(1,9)=1;
% Inductor
MNA(8,2)=1;
MNA(8,8)=-Lmatch/delt;
MNA(2,8)=1;
end

function [MNA,RHS,LHS,LHSpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current,vd,I]=MySPICE(MNA,RHS,LHS,LHSpre,Cdiode,Cload,Rload,Cmatch,Lmatch,Csrc,V0,V1,V2,G0,I0,tol_value,input_source,time,delt,delt_old,index,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre)
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
    Vnode10pre=LHS(1)-0;
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
            (-Lmatch/delt)*ILmatchpre;...
            (Csrc/delt)*Vnode10pre;];

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
%Pin=abs(fft_Vin.^2)/50;
%figure(3),hold on;plot(ff/1e9,(abs(Pin)));grid on;xlim([-30,30])

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
pout_val=abs(interp1(ff,(Pout),0,'nearest'));

% efficiency calculation
eff=(pout_val/pin_val)*100;
disp(['RF-DC conversion efficiency: ',num2str(eff),'%']);
end


function [MNA,RHS,LHS,Rload,Rser,Cload,Csrc]=SPICE_INIT_Rser_RC_NW(delt,delx,dely,delz)

MNA=zeros(4);

% load
Cload=1e-12;
Rload=25;
Rser=25;
cap_area=dely*delx;
cap_height=delz*2;
Csrc=8.854e-12*4.4*(cap_area)/(cap_height);

MNA(1,1)=1/Rser;
MNA(2,2)=1/Rser;
MNA(1,2)=-1/Rser;
MNA(2,1)=-1/Rser;

MNA(4,1)=Csrc/delt;
MNA(4,4)=-1;
MNA(1,4)=1;

% MNA matrix for the load resistor
MNA(2,2)=MNA(2,2)+1/Rload;

% MNA matrix for the load section
MNA(3,2)=Cload/delt;
MNA(3,3)=-1;
MNA(2,3)=1;

RHS=[0;0;0;0;];
LHS=[0;0;0;0;];

end

function [MNA,RHS,LHS]=MySPICE2(MNA,RHS,LHS,V1,V2,SRC,delt,delt_old,Rload,Rser,Cload,Csrc,time)
%  Vnode20pre=LHS(2);
%     Vnode10pre=LHS(1);
for n=time:delt:time+(delt_old-delt)
   Vnode20pre=LHS(2);
    Vnode10pre=LHS(1);
    RHS=[SRC;...
        0;...
        (Cload/delt)*Vnode20pre;...
        (Csrc/delt)*Vnode10pre;];

    LHS=MNA\RHS;
end
end

