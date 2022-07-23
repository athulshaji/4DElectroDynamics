%% Harmonic balance solver in time-domain
function MyHarmonicBalanceSolver_v2(freq)
% 5th Aug 2021
% Developed by Athul Shaji, Microwave Lab, ECE, IISc Bangalore
% LTSpice with HB for better results. Device models in LTSpice can be used
% in this way.
% In this code only 0 to K times the fundamental frequency are considered

% clc;close all;clear all;
format long;
clc;close all;
global qq;
global raw_data;
global frequencies freq_index K;
global v3 v2
global t
global vsrc;
global N;
global TIME START_TIME PAUSE_TIME counter
global f

dbm=5; % in dBm
SourceResValue=50; % Ohms
K=5;
inputpowerinwatt=10^(dbm/10)/1000;
Vrms=sqrt(inputpowerinwatt*SourceResValue);
Vpeak=Vrms*sqrt(2)*2;
TIME=100e-6+0.8163e-9;%0.8163e-9;%delt*61; % seconds
START_TIME=100e-6;
PAUSE_TIME=10;
omega0=freq; % fundamental freq
frequencies=[0:omega0:K*omega0];
% freq_index=[freq_index -1e9 -0.5e9 -0.25e9 0.25e9 0.5e9 1e9];
N=0;
% Fs=15*omega0;    
Fs=100*K*omega0; %1/(t(2)-t(1));
t=START_TIME:1/Fs:TIME; %linspace(0,TIME,512);
save t.txt t -ascii
raw_data=My_LTSpiceCall2(Vpeak,omega0,PAUSE_TIME);
var_mat=raw_data.variable_mat;
tt=raw_data.time_vect;
vn001=interp1(tt,var_mat(2,:),t,'spline');vsrc=vn001;clear VN001;
vn002=interp1(tt,var_mat(1,:),t,'spline');v2_pwl=[t.' real(vn002).' ];save v2_pwl.txt v2_pwl -ascii
vn003=interp1(tt,var_mat(3,:),t,'spline');v3_pwl=[t.' real(vn003).' ];save v3_pwl.txt v3_pwl -ascii
vsrc_pwl=[t.' vsrc.'];
save vsrc_pwl.txt vsrc_pwl -ascii
L=length(t);%L=L+1;
f = Fs*(0:L/2)/L;
% vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);
% f=Fs*(-L/2:L/2-1)/L;

VN002=abs((fft(vn002))/L);
VN002 = VN002(1:L/2+1);
VN002(2:end-1) = 2*VN002(2:end-1);

VN003=abs((fft(vn003))/L);
VN003 = VN003(1:L/2+1);
VN003(2:end-1) = 2*VN003(2:end-1);
% Finding the freq_indices
freq_index=[];
freq_index(1)=1;% DC
for p=1:length(frequencies)
   for q=1:length(f)
       if(abs((f(q)-frequencies(p))/frequencies(p))<1e-3)
           freq_index(p)=q;
       end
   end
end
a2=zeros(K+1+N,1);a2(1:K+1)=VN002(freq_index);
a3=zeros(K+1+N,1);a3(1:K+1)=VN003(freq_index);
clear v2 v3;

x0=[a2;a3]; % initial guess
options = optimoptions('fsolve',...
    'Algorithm','levenberg-marquardt',...
    'Display','iter-detailed',...    
    'FunctionTolerance',1e-3,...
    'StepTolerance',1e-3,...
    'MaxIterations',100000,...
    'MaxFunctionEvaluations',1000000,...
    'FiniteDifferenceType','central',...
    'ScaleProblem','jacobian',...
    'PlotFcn',[]);
% 
counter=1;
tic
[aoptim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
toc
exitflag
a2=aoptim(1:K+1+N);
a3=aoptim(K+1+N+1:end);
% alpha=aoptim(end-1:end);

for i=1:length(a2)
    v2(i,:)=a2(i).*exp(1i*2*pi.*((frequencies(i))).*(t));
end
for i=1:length(a3)
    v3(i,:)=a3(i).*exp(1i*2*pi.*((frequencies(i))).*(t));
end
v2=real(sum(v2));
save v2.txt v2 -ascii
v3=real(sum(v3));
save v3.txt v3 -ascii
raw_data=My_LTSpiceCall(PAUSE_TIME);
var_mat=raw_data.variable_mat;
tt=raw_data.time_vect;
ic1=interp1(tt,var_mat(4,:),t,'spline');
id1=interp1(tt,var_mat(5,:),t,'spline');
ir2=interp1(tt,var_mat(6,:),t,'spline');
ir1=interp1(tt,var_mat(7,:),t,'spline');
save ir1.txt ir1 -ascii;
save id1.txt ir1 -ascii;
save ic1.txt ir1 -ascii;
save ir2.txt ir2 -ascii;
figure,plot(t,v2,'LineWidth',1);axis tight;grid on;
figure,plot(t,v3,'LineWidth',1);axis tight;grid on;
% figure,plot(t,v3,'LineWidth',1);axis tight;grid on;
end

%%
function F=CostFunction(a)
global freq_index frequencies;
global K;
global raw_data
global v3 v2
global t
global vsrc;
global N
global TIME START_TIME PAUSE_TIME counter
global f;

a2=a(1:K+1+N);
a3=a(K+1+N+1:end);
clear v2 v3
for i=1:K+1+N
    v2(i,:)=a2(i).*exp(1i*2*pi.*((frequencies(i))).*(t));
end
for i=1:K+1+N
    v3(i,:)=a3(i).*exp(1i*2*pi.*((frequencies(i))).*(t));
end
v2=real(sum(v2));v3=real(sum(v3));
v2_pwl=[t.' real(v2).' ];
v3_pwl=[t.' real(v3).' ];
save v2_pwl.txt v2_pwl -ascii
save v3_pwl.txt v3_pwl -ascii
raw_data=My_LTSpiceCall(PAUSE_TIME);
var_mat=raw_data.variable_mat;
tt=raw_data.time_vect;
ic1=interp1(tt,var_mat(4,:),t,'spline');
id1=interp1(tt,var_mat(5,:),t,'spline');
ir2=interp1(tt,var_mat(6,:),t,'spline');
ir1=interp1(tt,var_mat(7,:),t,'spline');
F=[abs(ir1+id1);...
    abs(id1-ic1+ir2);];
% save v2_pwl.txt v2_pwl -ascii
% save v3_pwl.txt v3_pwl -ascii
clear v2 v3 ic1 id1 ir1 ir2
end

%%
function raw_data = My_LTSpiceCall(PAUSE_TIME)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"
global TIME START_TIME
%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft6.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.asc\r\n'...
'R1 v2 N001 50\r\n'...
'R2 0 v3 1000\r\n'...
'C1 v3 0 150p\r\n'...
'V1 N001 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\vsrc_pwl.txt"\r\n'...
'V2 v2 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v2_pwl.txt"\r\n'...
'V3 v3 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v3_pwl.txt"\r\n'...
'D1 v2 v3 MyDiode\r\n'...
['.tran 0 ',num2str(TIME),' ',num2str(START_TIME),' \r\n']...
'.model MyDiode D(Is=5u Rs=20 N=1.05 tt=1e-11 Cjo=0.14p Vj=0.34 M=0.4 Fc=0.5 Bv=2 Eg=0.69)\r\n'...
'.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
'.backanno\r\n'...
'.end\r\n'];
% writing netlist to a file
fid = fopen(netlist,'w+');
fprintf(fid,code);
fclose(fid);

dos('LTSpiceCall.bat');
pause(PAUSE_TIME);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.raw');
end

function raw_data = My_LTSpiceCall2(Vpeak,omega0,PAUSE_TIME)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"
global TIME START_TIME
%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft7.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft7.asc\r\n'...
'R1 N002 N001 50\r\n'...
'R2 0 N003 1000\r\n'...
'C1 N003 0 150p\r\n'...
['V1 N001 0 SINE(0 ',num2str(Vpeak),' ',num2str(omega0), ' 0 0 90)\r\n']...
'D1 N002 N003 MyDiode\r\n'...
['.tran 0 ',num2str(TIME),' ',num2str(START_TIME),' \r\n']...
'.model MyDiode D(Is=5u Rs=20 N=1.05 tt=1e-11 Cjo=0.14p Vj=0.34 M=0.4 Fc=0.5 Bv=2 Eg=0.69)\r\n'...
'.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
'.backanno\r\n'...
'.end\r\n'];
% writing netlist to a file
fid = fopen(netlist,'w+');
fprintf(fid,code);
fclose(fid);

dos('LTSpiceCall2.bat');
pause(PAUSE_TIME);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft7.raw');
end