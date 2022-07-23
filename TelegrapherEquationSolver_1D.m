function TelegrapherEquationSolver_1D
% Telegrapher's equations are solved using the assumption of lossless tx
% line.
clc;clear all;close all;
format long
Nz=2001;
Tsteps=5000;
% parameters of the tx line
R=0;L=2500e-9;G=0;C=1e-9; % Lossless tx line
Lz=2e-2; % 0.05 m; length of the tx line

% parameters of the source and load
RL=50;
Z0=sqrt(L/C);
RS=Z0;

vp=1/sqrt(L*C);
S=1; % Courant factor
delz=Lz/((Nz)-1);
delt=S*delz/vp; % CFL criterion
% Initializing the voltage and current vectors
V=zeros(Nz,1);
I=zeros(Nz,1);

% Defining some important constants
beta_1=2*delt/(RS*C*delz);
beta_2=2*delt/(RL*C*delz);
beta_3=(2*delt)/(C*delz);
r=delt^2/(L*C*delz^2);
u=beta_3/2;
q=delt/(L*delz);
t=0:Tsteps-1;
freq_sine_wave=5.2e9;
%freq_gaussian_wave=1.2e9;
% vsrc=1*cos(2*pi*freq*t);
%tau=1/(2*freq_gaussian_wave);
%t0=3*tau;
% vsrc=50*exp(-((t-t0)/(tau)).^2).*sin(2*pi*freq_sine_wave*t);
vsrc=2*sin(2*pi*freq_sine_wave*t*delt);
% vsrc=1*exp(-((t-t0)/(tau)).^2);
IC=[0;0;0;0;];
%figure('Color','w')
tic
for n=1:Tsteps    
    % Voltage update
    % node 1
    %V(1)=(1-beta_1)*V(1)-beta_3*I(1)+beta_1*vsrc(n);
    V(1)=vsrc(n);
    % nodes from 2 to Nz-1
    k=2:Nz-1;
    V(k)=V(k)-u*(I(k)-I(k-1));
    
    %     V(Nz)=(1-beta_2)*V(Nz)+(beta_3)*I(Nz-1);

    raw_data_3=My_LTSpiceCall_2(I(Nz-1),delt,IC,L,C,delz);
    V(Nz)=raw_data_3.variable_mat(1,end);
    IC(1)=raw_data_3.variable_mat(1,end);
    IC(2)=raw_data_3.variable_mat(2,end);
    IC(3)=raw_data_3.variable_mat(3,end);
    IC(4)=raw_data_3.variable_mat(7,end);
    Vleft(n)=raw_data_3.variable_mat(2,end);
    Vright(n)=raw_data_3.variable_mat(3,end);
    % Current update
    k=1:Nz-1;
    I(k)=I(k)-q*(V(k+1)-V(k));
    
    %plot(V);ylim([-2,2]);drawnow;
    
    if(mod(n,1000)==0)
        save('13march2022/case1.mat');
    end
end
toc
end


function raw_data = My_LTSpiceCall_2(value,delt,IC,L,C,delz)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

netlist = ['D:\Research\LTspice\Draft9.net'];
code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
    ['I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    [';V1 N001 0 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    [';L1 N001 N002 ',num2str(L*delz/2),'\r\n']...
    [';C1 N002 0 ',num2str(C*delz/2),'\r\n']...
    'C2 N001 N002 177.61f\r\n'...
    'L2 N002 0 3.41n Rser=1p\r\n'...
    'D1 N002 N003 MyDiode\r\n'...
    'C3 N003 0 1p\r\n'...
    'R1 N003 0 500\r\n'...
    ['.ic V(N002)=',num2str(IC(2)),' V(N003)=',num2str(IC(3)),' I(L2)=',num2str(IC(4)),'\r\n']...
    '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    '.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
    ['.tran ', num2str(delt),'\r\n']...
    '.backanno\r\n'...
    '.end\r\n'];

% writing netlist to a file
[fid, message] = fopen(netlist,'w+');
if fid < 0
    error('Failed to open myfile because: %s', message);
end
fprintf(fid,code);
fclose(fid);

% [status,cmdout]=dos('LTSpiceCall.bat');
% pause(0.5);
% [status,cmdout]=dos('LTSpiceEnd.bat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
working_folder = sprintf('D:\\Research\\LTspice');
path_of_exe='C:\"Program Files"\LTC\LTspiceXVII\XVIIx64.exe';
system(sprintf('cd %s && %s -Run -b %s -j2 && cd ..',working_folder,path_of_exe,netlist));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_data = LTspice2Matlab('D:\\Research\\LTspice\\Draft9.raw');
end






function raw_data = My_LTSpiceCall_0(value,delt,IC)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft8.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft8.asc\r\n'...
    [';I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    ['V1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    ';C1 N001 0 1n\r\n'...
    ';R1 N002 0 100\r\n'...
    ';C2 N002 0 50p\r\n'...
    'D1 N001 N002 MyDiode\r\n'... .
    ['.ic V(N001)=',num2str(IC(1)),'V V(N002)=',num2str(IC(2)),'V\r\n']...
    [';.ic V(N001)=',num2str(IC(1)),'\r\n']...
    ';.model D D\r\n'...
    '.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
    '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    ['.tran ', num2str(delt),'\r\n']...
    '.backanno\r\n'...
    '.end\r\n'];
% code  = ['* C:\\Users\\MWLAB\\OneDrive - Indian Institute of Science\\Desktop\\LTSpiceCkts\\RC_circuit\\RC_circuit.asc\r\n'...
% 'I1 0 N001 {CURRENT_VALUE}\r\n'...
% 'C1 N001 0 31.166e-15\r\n'...
% 'R1 N002 N001 25\r\n'...
% 'R2 N002 0 25\r\n'...
% '.params CURRENT_VALUE=',num2str(current_val),'\r\n'...
% '.tran 1ns\r\n'...
% '.backanno\r\n'...
% '.end\r\n'];
%

% writing netlist to a file
fid = fopen(netlist,'w+');
fprintf(fid,code);
pause(0.5)
fclose(fid);

dos('LTSpiceCall.bat');
pause(0.5);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft8.raw');
end

function raw_data = My_LTSpiceCall_1(value,delt,IC)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft8.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft8.asc\r\n'...
    [';I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    ['V1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    ';C1 N001 0 1n\r\n'...
    'R1 N001 0 50\r\n'...
    ';C2 N002 0 50p\r\n'...
    ';D1 N001 N002 1N914\r\n'... .
    ['.ic V(N001)=',num2str(IC(1)),'V\r\n']...
    ';.model D D\r\n'...
    ';.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
    ';.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    ['.tran ', num2str(delt),'\r\n']...
    '.backanno\r\n'...
    '.end\r\n'];
% code  = ['* C:\\Users\\MWLAB\\OneDrive - Indian Institute of Science\\Desktop\\LTSpiceCkts\\RC_circuit\\RC_circuit.asc\r\n'...
% 'I1 0 N001 {CURRENT_VALUE}\r\n'...
% 'C1 N001 0 31.166e-15\r\n'...
% 'R1 N002 N001 25\r\n'...
% 'R2 N002 0 25\r\n'...
% '.params CURRENT_VALUE=',num2str(current_val),'\r\n'...
% '.tran 1ns\r\n'...
% '.backanno\r\n'...
% '.end\r\n'];
%

% writing netlist to a file
fid = fopen(netlist,'w+');
fprintf(fid,code);
pause(0.5)
fclose(fid);

dos('LTSpiceCall.bat');
pause(0.5);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft8.raw');
end

% code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft8.asc\r\n'...
% ';I1 0 N001 ', num2str(current), '\r\n'...
% 'I1 0 N001 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\current_source.txt"\r\n'...
% 'C1 N001 0 3.3645f\r\n'...
% ';R1 N001 0 50\r\n'...
% 'D1 N001 0 MyDiode\r\n'...
% 'L1 N001 N002 47n\r\n'...
% 'C2 N002 0 100n\r\n'...
% 'R1 N002 0 500\r\n'...
% '.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
% '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
% ['.tran ', num2str(delt),'\r\n']...
% '.backanno\r\n'...
% '.end\r\n'];