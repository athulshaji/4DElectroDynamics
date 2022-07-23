function CNFDTD_1D_SPICE
%Code for solving Maxwell's equations in 1D using CNFDTD. Combining with SPICE.
% 02March2022

clc;close all;clear;format short

S=5;
num_time_steps=1000;
vamp=0.6325;
Nz=2001;
Lz=2e-2;
freq=5.2e9;
c0=299792458;
delz=Lz/(Nz-1);
dely=1.6e-3;
epsr=4.4;mur=1;
epsz=epsr*ones(Nz,1)*8.854e-12;
muz=mur*4*pi*1e-7*ones(Nz,1);
sigmaz=0*ones(Nz,1);
delt_cfl=delz*epsr/c0;
delt=S*delt_cfl;
t=0:delt:(num_time_steps-1)*delt;
vsrc=vamp*sin(2*pi*freq*t);
fsrc=vsrc/dely;
famp=max(fsrc);
E=zeros(Nz,1);
H=zeros(Nz,1);

Dzb=(1/(2*delz))*spdiags([-ones(Nz,1),ones(Nz,1)],[-1,0],Nz,Nz);
Dzf=(1/(2*delz))*spdiags([-ones(Nz,1),ones(Nz,1)],[0,1],Nz,Nz);

F=[E;H];
LHS1=(-epsz/delt).*speye(Nz);
LHS2=Dzb;
LHS3=Dzf;
LHS4=(-muz/delt).*speye(Nz);

LHS=[LHS1 LHS2;LHS3 LHS4];
[L,U,P]=lu(LHS);

RHS1=(-epsz/delt).*speye(Nz);
RHS2=-Dzb;
RHS3=-Dzf;
RHS4=(-muz/delt).*speye(Nz);

RHS=[RHS1 RHS2;RHS3 RHS4];
IC=[0;0;0;0;];
field=zeros(2*Nz,1);
%figure(1)
tic
for n=1:num_time_steps 
    %field(1)=fsrc(n);
    F(1)=F(1)+S*fsrc(n);%2*S*fsrc(n)+F(1);
    F=(U\(L\(P\(RHS*F))));    
    Ivalue=F(2*Nz)*dely;
    raw_data_2 = My_LTSpiceCall_2(Ivalue,delt,IC);
    Vvalue=raw_data_2.variable_mat(1,end);
    F(Nz)=Vvalue/dely;
    IC(1)=(raw_data_2.variable_mat(1,:));% Voltage @ node1
    IC(2)=(raw_data_2.variable_mat(2,:));% Voltage @ node2
    IC(3)=(raw_data_2.variable_mat(3,:));% Voltage @ node3
    IC(4)=(raw_data_2.variable_mat(7,:));% I(L2)
    Idiode(n)=raw_data_2.variable_mat(6,end);    
    Vleft(n)=IC(2);
    Vright(n)=IC(3);
    Vin(n)=vsrc(n);
    Iin(n)=F(Nz+2)*dely;
%         plot(F(2:Nz)*dely);grid on;
%         title(['time step: ', num2str(n)]);
%         %ylim([-2*famp-0.1,2*famp+0.1]);
%         ylim([-2*vamp-0.1,2*vamp+0.1]);
%         drawnow;
if(mod(n,1000)==0)
    save 21March2022\case-S=2.5-0dBm-0.6325.mat
end
end
toc
My_MATLAB_Conversion_Efficiency(Vleft,Vright,Idiode,delt,freq);
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

function raw_data = My_LTSpiceCall_2(value,delt,IC)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

netlist = ['D:\Research\LTspice\Draft9.net'];
code = ['* D:\\Research\\LTspice\\Draft9.asc\r\n'...
    ['I1 0 N001 PULSE(0 ',num2str(value),' 0 0 0 ',num2str(delt),' ', num2str(delt),')\r\n']...    
    'C2 N001 N002 177.61f\r\n'...
    'L2 N002 0 3.41n Rser=1p\r\n'...
    'C3 N003 0 1p\r\n'...
    'R1 N003 0 500\r\n'...
    'D1 N002 N003 MyDiode\r\n'...
    ['.ic V(N001)=',num2str(IC(1)),' V(N002)=',num2str(IC(2)),'V',' V(N003)=',num2str(IC(3)),'V',...
    ' I(L2)=',num2str(IC(4)),'\r\n']...
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
system(sprintf('cd %s && %s -Run -b %s -j6 && cd ..',working_folder,path_of_exe,netlist));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_data = LTspice2Matlab('D:\\Research\\LTspice\\Draft9.raw');
end
