function TxLine_CNFDTD_ver2_12March2022
clc;close all;clear;format short;
Tsteps=2500;
S=1;
freq=5.2e9;
Nz=2001;
R=0;L=2500e-9;G=0;C=1e-9; % Lossless tx line
Rs=50;
Lz=2e-2; % 0.05 m; length of the tx line
delz=Lz/(Nz-1);
vp=1/sqrt(L*C);
delt_cfl=delz/vp; % CFL criterion
delt=delt_cfl*S;
simulation_time=Tsteps*delt;
disp(simulation_time);

AMP=2; %dbmtov(Pin,Zin);Pin=5dBm;Zin=50ohms
t=0:Tsteps-1;
VSRC(t+1)=AMP*0.5*sin(2*pi*freq*t*delt);
IC=[0;0;0;0;];
Mz=Nz-2;



Dzb=(1/(2*delz))*spdiags([-ones(Mz,1),ones(Mz,1)],[-1,0],Mz,Mz);
Dzf=(1/(2*delz))*spdiags([-ones(Mz,1),ones(Mz,1)],[0,1],Mz,Mz);

LHS_1=(C/delt)*speye(Mz);
LHS_2=Dzb;
LHS_3=Dzf;
LHS_4=(L/delt)*speye(Mz);

RHS_1=(C/delt)*speye(Mz);
RHS_2=-Dzb;
RHS_3=-Dzf;
RHS_4=(L/delt)*speye(Mz);

LHS=[LHS_1 LHS_2;LHS_3 LHS_4];
tic;[L_,U_,P_]=lu(LHS);toc;
RHS=[RHS_1 RHS_2;RHS_3 RHS_4];
% Ipre1=zeros(Tsteps,1);Ipre2=zeros(Tsteps,1);

RHS_additional=[zeros(Mz,1);zeros(Mz,1)];
VI=[zeros(Mz,1);zeros(Mz,1)];

IC=[0;0;0;0;];
V2pre=0;Ickt=0;V1pre=0;V1=0;
tic
for n=1:Tsteps
    %I1=(2*delz)*((C/delt)*V1pre-(1/(2*delz))*I1pre-(C/delt)*VSRC(n));
    V1=VSRC(n);
    %V1=((2*S*delt)/(Rs*C*delz))*VSRC(n);
    RHS_additional(1)=S*(C/delt)*(V1pre-V1);%(I1+I1pre)/(2*delz);%I1/(2*delz);S*(C/delt)*VSRC(n)
    %V1pre=V1;
     VI=(U_\(L_\(P_*(RHS*VI+RHS_additional))));
    raw_data_3=My_LTSpiceCall_3((VI(2*Mz)+Ickt)/2,delt,IC,L,C,delz);
    %raw_data_3=My_LTSpiceCall_2((VI(2*Mz)+Ickt)/2,delt,IC,L,C,delz);
    IC(1)=raw_data_3.variable_mat(1,end);
    IC(2)=raw_data_3.variable_mat(2,end);
%     IC(3)=raw_data_3.variable_mat(3,end);
%     IC(4)=raw_data_3.variable_mat(7,end);
    RHS_additional(2*Mz)=-((raw_data_3.variable_mat(1,end)+V2pre))/(2*delz);    
    Ickt=raw_data_3.variable_mat(3,end);
    V2pre=raw_data_3.variable_mat(1,end);    
    
    Vnode1(n)=IC(1);
    Vnode2(n)=IC(2);

    if(mod(n,Tsteps)==0)
        save 12Apr2022\case-LTspice-S=1-timesteps=2500-Nz=2001-series-Res=50ohms-Cap=1pF.mat
    end

end
toc
My_MATLAB_Conversion_Efficiency(Vleft,Vright,Idiode,delt,freq);
end


function raw_data = My_LTSpiceCall_2_DC_Analysis(value)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

netlist = ['D:\Research\LTspice\Draft9.net'];
code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
    ['V1 N001 0 ',num2str(value),'\r\n']...
    'C2 N001 N002 177.61f\r\n'...
    'L2 N002 0 3.41n Rser=1p\r\n'...
    'D1 N002 N003 MyDiode\r\n'...
    'C3 N003 0 1p\r\n'...
    'R1 N003 0 500\r\n'...
    '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    '.lib C:\\Users\\MWLAB\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
    '.op\r\n'...
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
system(sprintf('cd %s && %s -Run -b %s -j6 && cd ..',working_folder,path_of_exe,netlist));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_data = LTspice2Matlab('D:\\Research\\LTspice\\Draft9.raw');
end

function raw_data = My_LTSpiceCall_2(value,delt,IC,L,C,delz)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

netlist = ['D:\Research\LTspice\Draft9.net'];
code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
    ['I1 0 N001 PULSE(',num2str(0),' ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    [';V1 N001 0 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
    [';L1 N001 N002 ',num2str(L*delz/2),'\r\n']...
    [';C1 N002 0 ',num2str(C*delz/2),'\r\n']...
    'C2 N001 N002 156.373f\r\n'...
    'L2 N002 0 3.83031n Rser=1p\r\n'...
    'D1 N002 N003 MyDiode\r\n'...
    'C3 N003 0 1p\r\n'...
    'R1 N003 0 832.98\r\n'...
    ['.ic V(N002)=',num2str(IC(2)),' V(N003)=',num2str(IC(3)),' I(L2)=',num2str(IC(4)),'\r\n']...
    '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    '.lib C:\\Users\\MWLAB\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
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
system(sprintf('cd %s && %s -Run -b %s -j6 && cd ..',working_folder,path_of_exe,netlist));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_data = LTspice2Matlab('D:\\Research\\LTspice\\Draft9.raw');
end

function raw_data=My_LTSpiceCall_3(value,delt,IC,L,C,delz)
netlist = ['D:\Research\LTspice\Draft9.net'];
code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
    ['I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...  
    'C1 N001 0 1n\r\n'...
    'R1 N001 N002 50\r\n'...
    ';D1 N001 0 MyDiode\r\n'...
    ';.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
    ';.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
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
system(sprintf('cd %s && %s -Run -b %s -j6 && cd ..',working_folder,path_of_exe,netlist));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_data = LTspice2Matlab('D:\\Research\\LTspice\\Draft9.raw');
end

function My_MATLAB_Conversion_Efficiency(Vleft,Vright,Idiode,delt,freq)
Vleft=Vleft(1:end);
Vright=Vright(1:end);
Idiode=Idiode(1:end);
NN=1e7;
Vleft=[Vleft,zeros(1,NN)];
Idiode=[Idiode,zeros(1,NN)];
Vright=[Vright,zeros(1,NN)];
LL=length(Vleft);
START=1;
STOP=length(Vleft);
INPUT=Vleft(START:STOP); %(V_antenna_port/dely).^2/(2*120*pi);
% Input
fft_Vin=fftshift(fft(INPUT))/LL;
%fft_Vin=fftshift(fft(V_txline_sample));
LL=length(fft_Vin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(1),hold on;plot(ff/1e9,log10(abs(fft_Vin)));grid on;xlim([-10,10])
fft_Iin=fftshift(fft(Idiode(START:STOP)))/LL;
LL=length(fft_Iin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(2),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pin=0.5*real(fft_Vin .* conj(fft_Iin));
%Pin=abs(fft_Vin.^2)/50;
 %figure(3),hold on;plot(ff/1e9,(abs(Pin)));grid on;xlim([-30,30])

% Output
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


% eff calc
pin_val=abs(interp1(ff,(Pin),freq*0,'nearest')) +...
    2*(abs(interp1(ff,(Pin),freq*1,'nearest')) +...
    abs(interp1(ff,(Pin),freq*2,'nearest')) +...
    abs(interp1(ff,(Pin),freq*3,'nearest')) +...
    abs(interp1(ff,(Pin),freq*4,'nearest')));
%pin_val=2*(abs(interp1(ff,(Pin),2.4e9,'nearest')));
pout_val=abs(interp1(ff,(Pout),0,'nearest'));

%pin_val=(10^(-10/10))/1000;

eff=(pout_val/pin_val)*100;
disp(eff);

end

% code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
%     ['I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
%     [';V1 N001 0 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...
%     'L1 N001 N002 12.5p\r\n'...
%     'C1 N002 0 5f\r\n'...
%     'C2 N002 N003 177.61f\r\n'...
%     'L2 N003 0 3.41n Rser=1p\r\n'...
%     'C3 N004 0 1p\r\n'...
%     'R1 N004 0 500\r\n'...
%     'D1 N003 N004 MyDiode\r\n'...
%     ['.ic V(N001)=',num2str(IC(1)),'V',' V(N002)=',num2str(IC(2)),'V',' V(N003)=',num2str(IC(3)),'V',' V(N004)=',num2str(4),'V',' I(L1)=',num2str(IC(5)),' I(L2)=',num2str(IC(6)),'\r\n']...
%     '.model MyDiode D(Is=5u Rs=20 N=1.05 TT=1e-11 Cjo=0.14p Vj=0.34 M=0.4 FC=0.5 Bv=2 IBV=1e-4 Xti=2 Eg=0.69)\r\n'...
%     '.lib C:\\Users\\MWLAB\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
%     ['.tran ', num2str(delt),'\r\n']...
%     '.backanno\r\n'...
%     '.end\r\n'];
%
% [';.ic V(N002)=',num2str(IC(1)),'V',' V(N003)=',num2str(IC(2)),'V',' V(N004)=',num2str(IC(3)),'V',' I(L1)=',num2str(IC(4)),...
%     ' I(L2)=',num2str(IC(5)),' I(L3)=',num2str(IC(6)),'\r\n']...

%
% VI(1)=VSRC(n);
%     raw_data_2=My_LTSpiceCall_2(VI((2*Nz)-1),delt,IC);
%     VI(Nz)=raw_data_2.variable_mat(2,end);
%
%     IC(1)=raw_data_2.variable_mat(2,end);% Voltage @ node2
%     IC(2)=raw_data_2.variable_mat(3,end);% Voltage @ node3
%     IC(3)=raw_data_2.variable_mat(4,end);% Voltage @ node4
%     IC(4)=raw_data_2.variable_mat(10,end);% I(L1)
%     IC(5)=raw_data_2.variable_mat(9,end);% I(L2)
%     IC(6)=raw_data_2.variable_mat(8,end);% I(L3)
%
%     Idiode(n)=raw_data_2.variable_mat(7,end);
%     Ires(n)=raw_data_2.variable_mat(12,end);
%     Icap(n)=raw_data_2.variable_mat(5,end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LHS_1(1,1)=1;LHS_2(1,1)=0;
% LHS_1(Nz,Nz)=1;%LHS_2(Nz,Nz)=0;%LHS_2(Nz,Nz-1)=0;
% %LHS_3(Nz,Nz)=0;
% %LHS_4(Nz,Nz)=1;
%
%
% LHS=[LHS_1 LHS_2;LHS_3 LHS_4];
% tic;[L_,U_,P_]=lu(LHS);toc;
%
% RHS_1=(C/delt)*speye(Nz);
% RHS_2=-Dzb;
% RHS_3=-Dzf;
% RHS_4=(L/delt)*speye(Nz);
%
% RHS_1(1,1)=1;RHS_2(1,1)=0;
% RHS_1(Nz,Nz)=1;%RHS_2(Nz,Nz)=0;%RHS_2(Nz,Nz-1)=0;
% %RHS_3(Nz,Nz)=0;
% %RHS_4(Nz,Nz)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[';.ic V(N001)=',num2str(IC(1)),'V',' V(N002)=',num2str(IC(2)),'V',' V(N003)=',num2str(IC(3)),'V',' I(L2)=',num2str(IC(4)),'\r\n']...
