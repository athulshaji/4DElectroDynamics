%% Harmonic balance solver in time-domain
function MyHarmonicBalanceSolver(freq)
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
global TIME

dbm=0; % in dBm
SourceResValue=50; % Ohms
K=2*5;
inputpowerinwatt=10^(dbm/10)/1000;
Vrms=sqrt(inputpowerinwatt*SourceResValue);
Vpeak=Vrms*sqrt(2)*2;
TIME=0.8163e-9;%delt*61; % seconds
omega0=freq; % fundamental freq
frequencies=[0:omega0:K*omega0];
% freq_index=[freq_index -1e9 -0.5e9 -0.25e9 0.25e9 0.5e9 1e9];
N=0;
% Fs=15*omega0;    
t=linspace(0,TIME,1024);
save t.txt t -ascii
Fs=1/(t(2)-t(1));
L=length(t);
f = Fs*(0:(L/2))/L;
freq_index=zeros(K+1,1);
freq_index(1)=1;% index corresponding to DC
for p=2:K+1
    for m=2:length(f)
        if(abs((f(m)-frequencies(p))/frequencies(p))<1e-3)
            freq_index(p)=m;
        end
    end
end

vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);

vsrc_pwl=[t.' vsrc];
save vsrc_pwl.txt vsrc_pwl -ascii
a2=zeros(K+1+N,1);
a3=zeros(K+1+N,1);clear v2 v3;
% a2(K)=0.588/2;
% a2(K+2)=0.588/2;
% a3(K+1)=0.465;
% for i=1:K+1+N
%     v2(i,:)=a2(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
% end
% for i=1:K+1+N
%     v3(i,:)=a3(i).*exp(1i*2*pi.*((freq_index(i))).*(t));
% end
x0=[a2;a3]; % initial guess
options = optimoptions('fsolve',...
    'Algorithm','levenberg-marquardt',...
    'Display','iter-detailed',...    
    'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4,...
    'MaxIterations',100000,...
    'MaxFunctionEvaluations',1000000,...
    'FiniteDifferenceType','central',...
    'PlotFcn',[]);
% 
tic
[aoptim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
toc
exitflag
a2=aoptim(1:K+1+N);
a3=aoptim(K+1+N+1:end);
% a3=aoptim(2*K+2+2*K+1:end);

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
raw_data=My_LTSpiceCall(0.7);
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
global TIME

a2=a(1:K+1+N);
a3=a(K+1+N+1:end);clear v2 v3;
% a2([1,3:end])=0;
% a3(2:end)=0;
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
raw_data=My_LTSpiceCall(0.7);
var_mat=raw_data.variable_mat;
tt=raw_data.time_vect;
ic1=interp1(tt,var_mat(4,:),t,'spline');
id1=interp1(tt,var_mat(5,:),t,'spline');
ir2=interp1(tt,var_mat(6,:),t,'spline');
ir1=interp1(tt,var_mat(7,:),t,'spline');
% ic1=interp1(tt,var_mat(8,:),t,'spline');
% id1=interp1(tt,var_mat(9,:),t,'spline');
% ir2=interp1(tt,var_mat(10,:),t,'spline');
% ir1=interp1(tt,var_mat(11,:),t,'spline');
% Fs=1/(t(2)-t(1));
% L=length(t);
% FREQ=Fs*(-L/2:L/2-1)/L;
% for m=1:K+1+N
%     v2_pwl=[t.' real(v2(m,:)).' ];
%     v3_pwl=[t.' real(v3(m,:)).' ];
%     save v2_pwl.txt v2_pwl -ascii
%     save v3_pwl.txt v3_pwl -ascii
%     raw_data=My_LTSpiceCall(0.262);
%     var_mat=raw_data.variable_mat;
%     tt=raw_data.time_vect;
%     ic1(m,:)=interp1(tt,var_mat(8,:),t);
%     id1(m,:)=interp1(tt,var_mat(9,:),t);
%     ir2(m,:)=interp1(tt,var_mat(10,:),t);
%     ir1(m,:)=interp1(tt,var_mat(11,:),t);
% end
% Y=fft(ir1+id1);% KCL for NODE1
% L=length(Y);
% P2 = (Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% NODE1=P1(freq_index).';
% % Fs=1/(t(2)-t(1));
% % f = Fs*(0:(L/2))/L;
% % stem(f,P1) 
% Y=fft(ir2-id1+ic1);% KCL for NODE2
% L=length(Y);
% P2 = (Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% NODE2=P1(freq_index).';
F=[(ir1+id1);...
    (id1+(ic1+ir2));];
% F=[NODE1;...
%     NODE2;];
CURRENT_NODE1=(ir1+id1);
save CURRENT_NODE1.txt CURRENT_NODE1 -ascii
CURRENT_NODE2=(ir2+id1+ic1);
save CURRENT_NODE2.txt CURRENT_NODE2 -ascii
clear v2 v3 ic1 id1 ir1 ir2
end

%%
function raw_data = My_LTSpiceCall(ic_v3)
%% USEFUL PATHS

% C:\Users\athul\OneDrive - Indian Institute of
% Science\Dev\MATLAB_Code\MATLAB_FDTD_SPICE\Diode_ckt.net

% "C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"
global TIME
%% CODE
netlist = ['D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\LTSpice\Draft6.net'];
code = ['* D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.asc\r\n'...
'R1 v2 N001 50\r\n'...
'R2 0 v3 1000\r\n'...
'C1 v3 0 150p\r\n'...
';V1 N001 0 SINE(0 0.632 2G 0 0 90)\r\n'...
'V1 N001 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\vsrc_pwl.txt"\r\n'...
'V2 v2 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v2_pwl.txt"\r\n'...
'V3 v3 0 PWL file="D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\v3_pwl.txt"\r\n'...
'D1 v2 v3 MyDiode\r\n'...
';D1 N002 N003 MyDiode\r\n'...
';B1 N002 0 V=V(v2)\r\n'...
';B2 N003 0 V=V(v3)\r\n'...
[';.ic V(C1)=',num2str(ic_v3),'\r\n']...
';.ic V(C1)=0.1\r\n'...
['.tran 0 ',num2str(TIME),' \r\n']...
'.model MyDiode D(Is=5u Rs=20 N=1.05 tt=1e-11 Cjo=0.14p Vj=0.34 M=0.4 Fc=0.5 Bv=2 Eg=0.69)\r\n'...
'.lib D:\\OneDrive - Indian Institute of Science\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio\r\n'...
'.backanno\r\n'...
'.end\r\n'];
% writing netlist to a file
fid = fopen(netlist,'w+');
fprintf(fid,code);
fclose(fid);

dos('LTSpiceCall.bat');
pause(1);
dos('LTSpiceEnd.bat');
raw_data = LTspice2Matlab('D:\\OneDrive - Indian Institute of Science\\Dev\\MATLAB_Code\\MATLAB_Code-master\\new_codes\\LTSpice\\Draft6.raw');
end