function [   tdjyxp, tdjzxp, tdmyxp, tdmzxp,...
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
    fdjxzn, fdjyzn, fdmxzn, fdmyzn]     = init_TD_FD_J_M(   freqs_of_interest,...
                                                            ffield_is, ffield_js, ffield_ks, ...
                                                            ffield_ie, ffield_je, ffield_ke)

% For +x face we have Jy,Jz, My,Mz according to J = x x H ; and M = -x x E;
tdjyxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current y component array calculated at +x surface
tdjzxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at +x surface
tdmyxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current y component array calculated at +x surface
tdmzxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at +x surface

% For -x face we have Jy,Jz, My,Mz according to J = -x x H ; and M = -(-x) x E;
tdjyxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current y component array calculated at -x surface
tdjzxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at -x surface
tdmyxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current y component array calculated at -x surface
tdmzxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at -x surface

%------------------------------------------------------------------------------------------------------------------------------------------

% For +y face we have Jx,Jz,Mx,Mz according to J = y x H ; and M = -y x E;
tdjxyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current x component array calculated at +y surface
tdjzyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at +y surface
tdmxyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current x component array calculated at +y surface
tdmzyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at +y surface

% For -y face we have Jx,Jz,Mx,Mz according to J = -y x H ; and M = -(-y) x E;
tdjxyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current x component array calculated at -y surface
tdjzyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain electric current z component array calculated at -y surface
tdmxyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current x component array calculated at -y surface
tdmzyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Time domain magnetic current z component array calculated at -y surface

%------------------------------------------------------------------------------------------------------------------------------------------

% For +z face we have Jx,Jy, Mx,My according to J = z x H ; and M = -z x E;
tdjxzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current x component array calculated at +z surface
tdjyzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current y component array calculated at +z surface
tdmxzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current x component array calculated at +z surface
tdmyzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current y component array calculated at +z surface

% For -z face we have Jx,Jy, Mx,My according to J = -z x H ; and M = -(-z) x E;
tdjxzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current x component array calculated at -z surface
tdjyzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain electric current y component array calculated at -z surface
tdmxzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current x component array calculated at -z surface
tdmyzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Time domain magnetic current y component array calculated at -z surface

%------------------------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------End of time domain electric andmagnetic current arrays-----------------------------------------

%--------------------------------------------Frequency domain electric and magnetic current arrays-----------------------------------------

% The first dimension in these array reserved for DFT frequencies

%------------------------------------------------------------------------------------------------------------------------------------------

% For +x face we have Jy,Jz, My,Mz according to J = x x H ; and M = -x x E;
fdjyxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current y component array calculated at +x surface
fdjzxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at +x surface
fdmyxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current y component array calculated at +x surface
fdmzxp=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at +x surface

% For -x face we have Jy,Jz, My,Mz according to J = -x x H ; and M = -(-x) x E;
fdjyxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current y component array calculated at -x surface
fdjzxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at -x surface
fdmyxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current y component array calculated at -x surface
fdmzxn=zeros(size(freqs_of_interest,2),1,ffield_je-ffield_js,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at -x surface

%----------------------------------------------------------------------------------------------------------------------------------------

% For +y face we have Jx,Jz,Mx,Mz according to J = y x H ; and M = -y x E;
fdjxyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current x component array calculated at +y surface
fdjzyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at +y surface
fdmxyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current x component array calculated at +y surface
fdmzyp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at +y surface

% For -y face we have Jx,Jz,Mx,Mz according to J = -y x H ; and M = -(-y) x E;
fdjxyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current x component array calculated at -y surface
fdjzyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain electric current z component array calculated at -y surface
fdmxyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current x component array calculated at -y surface
fdmzyn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,1,ffield_ke-ffield_ks); % Frequency domain magnetic current z component array calculated at -y surface

%--------------------------------------------------------------------------------------------------------------------------------------

% For +z face we have Jx,Jy, Mx,My according to J = z x H ; and M = -z x E;
fdjxzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current x component array calculated at +z surface
fdjyzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current y component array calculated at +z surface
fdmxzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current x component array calculated at +z surface
fdmyzp=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current y component array calculated at +z surface

% For -z face we have Jx,Jy, Mx,My according to J = -z x H ; and M = -(-z) x E;
fdjxzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current x component array calculated at -z surface
fdjyzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain electric current y component array calculated at -z surface
fdmxzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current x component array calculated at -z surface
fdmyzn=zeros(size(freqs_of_interest,2),ffield_ie-ffield_is,ffield_je-ffield_js,1); % Frequency domain magnetic current y component array calculated at -z surface

%--------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------End of frequency domain electric and magnetic current arrays--------------------------------

%-------------------------------------End of Near Field Far Field - Frequency Domain Transformation Arrays-----------------------------

end