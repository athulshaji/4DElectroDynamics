function CN_FDTD_Solver2(T, ...
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

FREQ = fmax;
VIDEO = 1; % numerals 0 or 1
PML_ENABLED = 1;
file_name = '21March2022.avi';
%% Basic parameter

c0 = 299792458; %m/s
Nx = nxtot;
Ny = nytot;
Nz = nztot;

mu0 = 4*pi*1e-7; % A/m
epsilon0 = 8.854e-12; % F/m
epsilonR = ones(Nx,Ny,Nz);
% epsilonR(ceil(Nx/2)-10:ceil(Nx/2)+5, ceil(Ny/2)-10:ceil(Ny/2)+5, ceil(Nz/2)-10:ceil(Nz/2)+5) = 1;
epsilon = (epsilon0.*epsilonR).*eps_x(1:Nx,1:Ny,1:Nz);
% sigma_x = zeros(Ny,Nx,Nz);
% sigma_y = zeros(Ny,Nx,Nz);
% sigma_z = zeros(Ny,Nx,Nz);

NxOriginal = Nx - 2*pml_thickness;
NyOriginal = Ny - 2*pml_thickness;
NzOriginal = Nz - 2*pml_thickness;

NxPML = pml_thickness;
NyPML = NxPML;
NzPML = NxPML;

% [NxOriginal, NyOriginal, NzOriginal]
% [Nx, Ny, Nz]

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

% [size_x,size_y,size_z] = size(Ex);
% Nxy = size_x*size_y;
% Nxyz = size_x*size_y*size_z;

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

sigma_x = sigma_x(1:Nx,1:Ny,1:Nz);


sigma_x_vec = [];
sigma_y_vec = [];
sigma_z_vec = [];
epsilon_vec = [];
for k = 1:Nz
    sigma_x_vec = [sigma_x_vec;reshape(sigma_x(:,:,k).',Nxy,1)];
    %     sigma_y_vec = [sigma_y_vec;reshape(sigma_y(:,:,k).',Nxy,1)];
    %     sigma_z_vec = [sigma_z_vec;reshape(sigma_z(:,:,k).',Nxy,1)];
    epsilon_vec = [epsilon_vec;reshape(epsilon(:,:,k).',Nxy,1)];
end

%% Defining the PML parameters
if(PML_ENABLED)
    m = 1; R = 1e-16;
    sigmaXOPT = -((m+1)*log(R)/(2*120*pi*NxPML*delx));
    sigmaXMAX = 20*sigmaXOPT;
    sigmaYOPT = -((m+1)*log(R)/(2*120*pi*NyPML*dely));
    sigmaYMAX = 20*sigmaYOPT;
    sigmaZOPT = -((m+1)*log(R)/(2*120*pi*NzPML*delz));
    sigmaZMAX = 20*sigmaZOPT;
    alphaXMax = 0.21; alphaYMax = alphaXMax; alphaZMax = alphaXMax;
    kappaXMax = 20; kappaYMax = kappaXMax; kappaZMax = kappaXMax;
    %sigmaXMax = 100; sigmaYMax = sigmaXMax; sigmaZMax = sigmaXMax;
    kappaX = ones(Ny,Nx,Nz); kappaY = kappaX; kappaZ = kappaX;
    
    for k = 1:Nz
        for j =1:Ny
            for i = 1:NxPML
                %alphaX(j,i,k) = alphaXMax * (i/NxPML)^m;
                alphaX(j,i,k) = alphaXMax * ((NxPML-i+1)/NxPML)^m;
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
                alphaY(j,i,k) = alphaYMax * ((NyPML - j + 1)/NyPML)^m;
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
                alphaZ(j,i,k) = alphaZMax * ((NzPML - k + 1)/NzPML)^m;
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
    %     CHx_current_vec = [];
    %     CHy_current_vec = [];
    %     CHz_current_vec = [];
    %     CHx_past_vec = [];
    %     CHy_past_vec = []
    %     CHz_past_vec = [];
    %     CEx_current_vec = []
    %     CEy_current_vec = [];
    %     CEz_current_vec = [];
    %     CEx_past_vec = [];
    %     CEy_past_vec = [];
    %     CEz_past_vec = [];
    ExIncVec = [];
    EyIncVec = [];
    EzIncVec = [];
    HxIncVec = [];
    HyIncVec = [];
    HzIncVec = [];
    
    for k = 1:Nz
        %         CHx_current_vec = [CHx_current_vec;reshape(CHx_current(:,:,k).',Nxy,1)];
        %         CHy_current_vec = [CHy_current_vec;reshape(CHy_current(:,:,k).',Nxy,1)];
        %         CHz_current_vec = [CHz_current_vec;reshape(CHz_current(:,:,k).',Nxy,1)];
        %
        %         CHx_past_vec = [CHx_past_vec;reshape(CHx_past(:,:,k).',Nxy,1)];
        %         CHy_past_vec = [CHy_past_vec;reshape(CHy_past(:,:,k).',Nxy,1)];
        %         CHz_past_vec = [CHz_past_vec;reshape(CHz_past(:,:,k).',Nxy,1)];
        %
        %         CEx_current_vec = [CEx_current_vec;reshape(CEx_current(:,:,k).',Nxy,1)];
        %         CEy_current_vec = [CEy_current_vec;reshape(CEy_current(:,:,k).',Nxy,1)];
        %         CEz_current_vec = [CEz_current_vec;reshape(CEz_current(:,:,k).',Nxy,1)];
        %
        %         CEx_past_vec = [CEx_past_vec;reshape(CEx_past(:,:,k).',Nxy,1)];
        %         CEy_past_vec = [CEy_past_vec;reshape(CEy_past(:,:,k).',Nxy,1)];
        %         CEz_past_vec = [CEz_past_vec;reshape(CEz_past(:,:,k).',Nxy,1)];
        
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
% a1 = dt/((2*epsilon));
% a2 = dt/(2*mu0);

[Dxe,Dxxe,Dxh,Dye,Dyye,Dyh,Dze,Dzze,Dzh] = DerivativeOperator(delx,dely,delz,Nx,Ny,Nz,Nxy,Nxyz);
a0 = (2.*epsilon_vec - sigma_x_vec.*dt)./(2*epsilon_vec + sigma_x_vec.*dt); % 1
a1 = dt./(2.*epsilon_vec);
a2 = dt./(2*mu0);
i = 1:Nxyz;
Dxew = ((Dxe(i,:).*a2.')./WX.');
% Dxxew = (Dxxe(i,:)./(WX.'));
Dxhw = ((Dxh(i,:).*a1.')./WX.');

Dyew = ((Dye(i,:).*a2.')./WY.');
% Dyyew = (Dyye(i,:)./(WY.'));
Dyhw = ((Dyh(i,:).*a1.')./WY.');

Dzew = ((Dze(i,:).*a2.')./WZ.');
% Dzzew = (Dzze(i,:)./(WZ.'));
Dzhw = ((Dzh(i,:).*a1.')./WZ.');
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
% A = speye(6*Nxyz,6*Nxyz);
A = spdiags([a0; a0; a0; ones(Nxyz,1); ones(Nxyz,1); ones(Nxyz,1)], 0,6*Nxyz, 6*Nxyz);
Z = sparse(Nxyz,Nxyz);

D1 = [  Z Z Z Z -Dzhw Z;...
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

PSI = [ RY.*psiEXY - RZ.*psiEXZ;...
    RZ.*psiEYZ - RX.*psiEYX;...
    RX.*psiEZX - RY.*psiEZY;...
    RZ.*psiHXZ - RY.*psiHXY;...
    RX.*psiHYX - RZ.*psiHYZ;...
    RY.*psiHZY - RX.*psiHZX;];

EYE = speye(6*Nxyz,6*Nxyz);

tic;
[L_STAR, U_STAR, P_STAR] = lu(EYE - D1);
[L,U,P] = lu(EYE - D2);
toc;

Z = sparse(Nxyz,1);
%% Output

if(VIDEO)
    vidObj = VideoWriter(file_name);
    open(vidObj);
end

INDEX = 0;
PEC_LOC = find(sigma_x == 1e20);

if(VIDEO)
    fig = figure('color','w','units','normalized','outerposition',[0 0 1 1]);
end

fieldEx=zeros(Nxyz,1);fieldEy=zeros(Nxyz,1);fieldEz=zeros(Nxyz,1);
fieldHx=zeros(Nxyz,1);fieldHy=zeros(Nxyz,1);fieldHz=zeros(Nxyz,1);
INPUT_FIELD=[fieldEx;fieldEy;fieldEz;fieldHx;fieldHy;fieldHz;];

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
    %             Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2);
    %             Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2);
    %             Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2);
    %             Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2);
    %             Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2);
    %             Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2);
    
    % Modulated Gaussian Source
%         Exincident_current = Exincamp * exp(-((((n-1)*delt+delt)-t0-krEx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
%         Eyincident_current = Eyincamp * exp(-((((n-1)*delt+delt)-t0-krEy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
%         Ezincident_current = Ezincamp * exp(-((((n-1)*delt+delt)-t0-krEz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt)-t0));
%         Hxincident_current = Hxincamp * exp(-((((n-1)*delt+delt/2)-t0-krHx)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
%         Hyincident_current = Hyincamp * exp(-((((n-1)*delt+delt/2)-t0-krHy)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
%         Hzincident_current = Hzincamp * exp(-((((n-1)*delt+delt/2)-t0-krHz)/tau).^2).*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0));
    
    
    % Sinusoidal Source
        Exincident_current = Exincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEx));
        Eyincident_current = Eyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEy));
        Ezincident_current = Ezincamp .*sin(2*pi*(FREQ)*(((n-1)*delt)-t0-krEz));
        Hxincident_current = Hxincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHx));
        Hyincident_current = Hyincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHy));
        Hzincident_current = Hzincamp .*sin(2*pi*(FREQ)*(((n-1)*delt+delt/2)-t0-krHz));
    
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
    %PHI_OLD = PHI_OLD+INC_FIELD ;
    RHS_STAR = (A+D1+(2*D2))*PHI_OLD + (PSI+INC_FIELD);
    PHI_STAR = U_STAR \ (L_STAR \ (P_STAR * RHS_STAR));
    RHS = PHI_STAR - (D2*(PHI_OLD));
    PHI_NEW = U \ (L \ (P * (RHS)));
    
    %PHI_NEW = PHI_NEW;
    
    ExNew = PHI_NEW(0*Nxyz+1:1*Nxyz,1);
    EyNew = PHI_NEW(1*Nxyz+1:2*Nxyz,1);
    EzNew = PHI_NEW(2*Nxyz+1:3*Nxyz,1);
    HxNew = PHI_NEW(3*Nxyz+1:4*Nxyz,1);
    HyNew = PHI_NEW(4*Nxyz+1:5*Nxyz,1);
    HzNew = PHI_NEW(5*Nxyz+1:6*Nxyz,1);

%     xloc=105;yloc=67;zloc=22;
%     xyzloc=Nxy*(zloc-1)+Ny*(xloc-1)+yloc;
%     fieldEx(xyzloc)=S*sin(2*pi*(FREQ)*(((n-1)*delt)-t0));
%     INPUT_FIELD=[fieldEx;fieldEy;fieldEz;fieldHx;fieldHy;fieldHz;];
    
    ExNew(PEC_LOC) = 0;
    EyNew(PEC_LOC) = 0;
    EzNew(PEC_LOC) = 0;
%     
%     ExNew = ExNew + ExIncVec;
%     EyNew = EyNew + EyIncVec;
%     EzNew = EzNew + EzIncVec;
    
    ExOld = PHI_OLD(0*Nxyz+1:1*Nxyz,1);
    EyOld = PHI_OLD(1*Nxyz+1:2*Nxyz,1);
    EzOld = PHI_OLD(2*Nxyz+1:3*Nxyz,1);
    HxOld = PHI_OLD(3*Nxyz+1:4*Nxyz,1);
    HyOld = PHI_OLD(4*Nxyz+1:5*Nxyz,1);
    HzOld = PHI_OLD(5*Nxyz+1:6*Nxyz,1);
    
%     ExOld(PEC_LOC) = 0;
%     EyOld(PEC_LOC) = 0;
%     EzOld(PEC_LOC) = 0;
       
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
    
    % vector field to cube
    Q = ExNew;
    for i = 1:Nz
        FIELDX(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Nx,Ny).');
%         FIELDX(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Ny,Nx));
    end
    
    Q = EyNew;
    for i = 1:Nz
        FIELDY(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Nx,Ny).');
        %         FIELDY(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Ny,Nx));
    end
    
    Q = EzNew;
    for i = 1:Nz
        FIELDZ(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Nx,Ny).');
        %         FIELDZ(:,:,i) = full(reshape(Q((i-1)*Nxy+1:i*Nxy,1),Ny,Nx));
    end
    
    %     FIELD_MAX(timeIter) = log10(max(QQ(:)));
    if(VIDEO)
        el = 180;
        az = 0;
        %         subplot 231
        %         slice(Exincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal
        %         subplot 232
        %         slice(Eyincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal
        %         subplot 233
        %         slice(Ezincident_current,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)),view(el,az),shading interp, title(['iter#',num2str(n)]);colorbar;axis equal
        
        val = 0.05;
        subplot 121 %234
        slice(Eyincident_current,ceil(Ny/2),ceil(Nx/2),ceil(Nz/2)+1),view(3); title(['iter#',num2str(n)]);
        caxis([-val, val]);
        colorbar;
        axis equal
        subplot 122 %235
        slice(FIELDY,ceil(Ny/2),ceil(Nx/2),ceil(Nz/2)+2),view(0,90); title(['iter#',num2str(n)]);
        caxis([-val, val]);
        colorbar;
        axis equal        
        %         subplot 133 %236
        %         slice(FIELDZ,ceil(Ny/2),ceil(Nx/2),ceil(Nz/2)),view(az,el),shading interp, title(['iter#',num2str(n)]);
        %         caxis([-val, val]);
        %colorbar;axis equal
        drawnow;
        FRAME = getframe(fig);
        writeVideo(vidObj, FRAME)
    end
end
toc
if(VIDEO)
    close(vidObj);
end
end