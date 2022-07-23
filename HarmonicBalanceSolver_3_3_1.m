function HarmonicBalanceSolver_3_3_1
% Developed by Athul Shaji, Microwave Lab, ECE, IISc Bangalore
% LTSpice with HB for better results. Device models in LTSpice can be used
% in this way.
%

clc;close all;clear all;
format long;

global qq;
global raw_data;
global freq_index K;
global v3 v2
global t
global vsrc;

dbm=0; % in dBm
SourceResValue=50; % Ohms
K=5;
inputpowerinwatt=10^(dbm/10)/1000;
Vrms=sqrt(inputpowerinwatt*SourceResValue);
Vpeak=Vrms*sqrt(2)*2;
T=1e-9; % seconds
omega0=2e9; % fundamental freq
freq_index=[-K*omega0:omega0:K*omega0];
Fs=100*omega0;    
t=linspace(0,T,Fs*T);
delt=t(20)-t(19);
vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);
% vsrc_pwl=[];
% for i=1:length(t)
%    vsrc_pwl=[vsrc_pwl t(i) vsrc(i)]; 
% end
vsrc_pwl=[t.' vsrc];
save vsrc_pwl.txt vsrc_pwl -ascii
% a2=[0;0;0;0;0.588847/2;0;0.588847/2;0;0;0;0];%rand(2*K+1,1);
% a3=[0;0;0;0;0;0.2;0;0;0;0;0];%rand(2*K+1,1);
a2=zeros(2*K+1,1);
a3=zeros(2*K+1,1);
% a2(K)=0.588847/2;
% a2(K+2)=0.588847/2;
% a3(K+1)=0.4632;
for i=1:2*K+1
    v2(i,:)=a2(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
for i=1:2*K+1
    v3(i,:)=a3(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
x0=[a2;a3]; % initial guess
options = optimoptions('fsolve',...
    'Display','iter-detailed',...
    'Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-8,...
    'StepTolerance',1e-8,...
    'MaxIterations',100000,...
    'MaxFunctionEvaluations',1000000,...
    'FiniteDifferenceType','central',...
    'PlotFcn',[]);
tic
[aoptim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
toc
exitflag
a2=aoptim(1:2*K+1);
a3=aoptim(2*K+2:end);
% a3=aoptim(2*K+2+2*K+1:end);

for i=1:2*K+1
    v2(i,:)=a2(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
for i=1:2*K+1
    v3(i,:)=a3(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
v2=real(sum(v2));
v3=real(sum(v3));
figure,plot(t,v2,'LineWidth',1);axis tight;grid on;
figure,plot(t,v3,'LineWidth',1);axis tight;grid on;
% figure,plot(t,v3,'LineWidth',1);axis tight;grid on;
end

function F=CostFunction(a)
global freq_index;
global K;
global raw_data
global v3 v2
global t
global vsrc;

% x0=a; % initial guess
% options = optimoptions('fsolve',...
%     'Display','iter-detailed',...
%     'Algorithm','levenberg-marquardt',...
%     'FunctionTolerance',1e-8,...
%     'StepTolerance',1e-8,...
%     'MaxIterations',100000,...
%     'MaxFunctionEvaluations',1000000,...
%     'FiniteDifferenceType','central',...
%     'PlotFcn',[]);
% [optim_value]=fsolve(@CostFunction2,x0,options);
% a=optim_value;
a2=a(1:2*K+1);
a3=a(2*K+2:end);
for i=1:2*K+1
    v2(i,:)=a2(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
for i=1:2*K+1
    v3(i,:)=a3(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
end
v2=real(sum(v2));v3=real(sum(v3));
v2_pwl=[t.' real(v2).' ];
v3_pwl=[t.' real(v3).' ];
fid=fopen('v2_pwl.txt');
fclose(fid);
fid=fopen('v3_pwl.txt');
fclose(fid);
save v2_pwl.txt v2_pwl -ascii
save v3_pwl.txt v3_pwl -ascii
raw_data=My_LTSpiceCall;
var_mat=raw_data.variable_mat;
tt=raw_data.time_vect;
ic1=interp1(tt,var_mat(4,:),t);
id1=interp1(tt,var_mat(5,:),t);
ir2=interp1(tt,var_mat(6,:),t);
ir1=interp1(tt,var_mat(7,:),t);
M=[id1+ir1;...
    (ir2+ic1)+id1;];
F=norm(M,'fro').^2;
clear v2 v3 ic1 id1 ir1 ir2
end

% function F=CostFunction2(x)
% global freq_index;
% global K;
% global raw_data
% global v3 v2
% global t
% global vsrc;
% 
% a=x;
% a2=a(1:2*K+1);
% a3=a(2*K+2:end);
% for i=1:2*K+1
%     v2(i,:)=a2(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
% end
% for i=1:2*K+1
%     v3(i,:)=a3(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
% end
% v2=real(sum(v2));v3=real(sum(v3));
% F=[v2(1)-v2(end);...
%     v3(1)-v3(end);];
% end

function raw_data = My_LTSpiceCall(freq)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"

%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft6.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.asc\r\n'...
'V1 N001 0 SINE(0 0.632 2e9 0 0 90)\r\n'...
'R1 N002 N001 50\r\n'...
'R2 0 N003 1000\r\n'...
'D1 N002 N003 MyDiode\r\n'...
'C1 N003 0 1u\r\n'...
'V2 N002 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v2_pwl.txt"\r\n'...
'V3 N003 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v3_pwl.txt"\r\n'...
'.tran 0 1n\r\n'...
'.model MyDiode D(Is=5u Rs=20 N=1.05 tt=1e-11 Cjo=0.14p Vj=0.34 M=0.4 Fc=0.5 Bv=2 Eg=0.69)\r\n'...
'.backanno\r\n'...
'.end\r\n'];
% 

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
fclose(fid);

dos('LTSpiceCall.bat');
pause(1);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.raw');
end

