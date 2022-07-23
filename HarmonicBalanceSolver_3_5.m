function HarmonicBalanceSolver_3_5
% Developed by Athul Shaji, Microwave Lab, ECE, IISc Bangalore
% Class A power amplifier. Class A --> Common emitter amplifier which can
% have a very large voltage gain
% Developing a time-domain harmonic balance
% i1a=linear current at node 1
% i1b=linear current at node 1
% i1_=nonlinear current at node 1
% i2=linear current at node 2
% i2_=nonlinear current at node 2
% i3=linear current at node 3
% i3_=nonlinear current at node 3
% WORKING!!!

clc;clear all;close all;

global lengthfreqvec;
global vsrc;
global positive_freq_index negative_freq_index freq_index;
global SeriesResistorValue CE RE;
global Fs;
global Rs;
global Isat;
global Vth;
global eta;
global K;
global f;
global SourceResValue
global Diode_Data
global REAL_DIODE_DATA DIODE_EQN_BASED
global Bv;
global omega0
global t
global vcc
global beta;
global R1 R2 C RC RL

% counter=1;
% fid=fopen('SMS7630_Schottky_Diode.txt');
% fgetl(fid);
% while(~feof(fid))
%     m=fgetl(fid);
%     m=regexp(m,'	','split');
%     m{1,1}=str2num(m{1,1});
%     m{1,2}=str2num(m{1,2});
%     Diode_Data.voltage(counter)=m{1};
%     Diode_Data.current(counter)=m{2};
%     counter=counter+1;
% end
% fclose(fid);
iter_number=1;
dbm_array=[-30];
tic
for dbm=dbm_array
    tic
    clc
    dbm
    vcc=5; % Volts
    beta=100;
    RL=500;
    R1=20e3;
    R2=3.6e3;
    RE=220;
    CE=1e-6;
    RC=1.2e3;
    C=0.01e-6;
    REAL_DIODE_DATA=0;
    DIODE_EQN_BASED=~REAL_DIODE_DATA;
    DO_PLOT=0;
    Bv=-2; %volts
    Rs=20; % Ohms
    Isat=5e-6;%2.52e-9; % Amperes for Si diode
    TempC=25; % degree Celcius
    TempK=TempC+273.15; % Kelvin
    Vth=1.38e-23*TempK/1.6e-19; %0.026 V for T=25degree Celcius;
    eta=1.05;%1.752;
    K=5; % number of harmonics of interest k=0:K
    SeriesResistorValue=50; % Ohms
    ShuntResValue=500; % Ohms
    ShuntCapValue=100e-9;
    inputpowerinwatt=10^(dbm/10)/1000;
    SourceResValue=50;
    % Irms=sqrt(inputpowerinwatt/SourceResValue);
    % Vrms=Irms*SourceResValue;
    Vrms=sqrt(inputpowerinwatt*SourceResValue);
    Vpeak=Vrms*sqrt(2)*2;
    T=0.04; % seconds
    omega0=50; % fundamental freq
    freq_index=[-K*omega0:omega0:K*omega0];
    Fs=100*omega0;    
    t=linspace(0,T,Fs*T);
    vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);
    der_vsrc=diff(vsrc)/(t(3)-t(2));
    der_vsrc=[vsrc(1);der_vsrc];
    L=length(vsrc);
    f=Fs*(-L/2:L/2-1)/L;
    % positive frequencies
    index=1;
    for i=0:omega0:K*omega0
        positive_freq_index(index,1)=find(f==i,1);
        index=index+1;
    end
    % negative frequencies
    index=1;
    for i=-K*omega0:omega0:-omega0
        negative_freq_index(index,1)=find(f==i,1);
        index=index+1;
    end
    freq_index=[];
    freq_index=[freq_index;negative_freq_index];
    freq_index=[freq_index;positive_freq_index];
    a1=zeros(2*K+1,1);
    a2=zeros(2*K+1,1);
    a3=zeros(2*K+1,1);
    global time K
    time=1:length(t);
    x0=[a1;a2;a3]; % initial guess
    options = optimoptions('fsolve',...
        'Display','iter-detailed',...
        'Algorithm','levenberg-marquardt',...
        'FunctionTolerance',1e-6,...
        'StepTolerance',1e-6,...
        'MaxIterations',100000,...
        'MaxFunctionEvaluations',1000000,...
        'FiniteDifferenceType','central',...
        'PlotFcn',[]);
    [a_optim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
    exitflag    
end
a1=a_optim(1:2*K+1);
a2=a_optim(2*K+2:2*K+2+2*K);
a3=a_optim(2*K+2+2*K+1:end);
for i=1:2*K+1
    v1(i,:)=a1(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));    
end
for i=1:2*K+1    
    v2(i,:)=a2(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));
end
for i=1:2*K+1    
    v3(i,:)=a3(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));
end
v1=real(sum(v1));v2=real(sum(v2));v3=real(sum(v3));
figure('Color','w');
subplot 311
plot(t,v1,'-ok','LineWidth',1);axis tight;grid on;
subplot 312
plot(t,v2,'-ok','LineWidth',1);axis tight;grid on;
subplot 313
plot(t,v3,'-ok','LineWidth',1);axis tight;grid on;
end

function F=CostFunction(a)
global time K
global vsrc;
global freq_index;
global CE RE;
global Isat;
global Vth;
global eta;
global K;
global f;
global Diode_Data
global Bv;
global omega0
global t
global vcc
global beta;
global R1 R2 RC RL C

a1=a(1:2*K+1);
a2=a(2*K+2:2*K+2+2*K);
a3=a(2*K+2+2*K+1:end);

for i=1:2*K+1
    v1(i,:)=a1(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));    
end
for i=1:2*K+1    
    v2(i,:)=a2(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));
end
for i=1:2*K+1    
    v3(i,:)=a3(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));
end
v1=(sum(v1));
v2=(sum(v2));
v3=(sum(v3));

% Node 1
i1a=(v1-vcc)/R1;
i1b=(v1-0)/R2;
i1c=diff(vsrc)/(t(3)-t(2));
i1c=C*[vsrc(1);i1c];
% i1c=C*dvdt(v1-vsrc);
vbe=(v1-v3);
for T=time
    if (vbe(T)>0)
        i1_(T)=Isat*(exp(vbe(T)/(eta*Vth))-1);
%         i1_(T)=interp1(real(Diode_Data.voltage),...
%             real(Diode_Data.current),real(vdiode(T)),'linear');
    elseif(vbe(T)<0&&vbe(T)>=Bv)
        i1_(T)=-Isat;
    elseif(vbe(T)<Bv)
        i1_(T)=interp1(real(Diode_Data.voltage),...
            real(Diode_Data.current),real(vbe(T)),'linear');
    elseif(vbe(T)==0)
        i1_(T)=0;
    end
end
i1total=i1a+i1b+i1c+i1_;
% Node 2
i2a=(v2-vcc)/RC;
i2b=(v2)/RL;
i2_= beta*i1_;% collector current
i2total=i2a+i2b+i2_;

% Node 3
i3_=(beta+1)*i1_; % emitter current
V3=fftshift(fft(v3))/length(v3);
V3=V3(freq_index);
Z=@(f) RE./(1+(1i*2*pi*f*RE*CE));
IMPEDANCE=Z((f(freq_index)));
I3=(V3)./IMPEDANCE;
I33=zeros(1,length(vsrc));
I33(freq_index)=I3;
i3=ifft(ifftshift((I33)))*length(vsrc);
i3total=i3+i3_;

% Cost function
F=[i1total;i2total;i3total];
clear v1 v2 v3
end

function derivative=dvdt(v,delt)
L=length(v);
for i=1:L-1
   derivative(i)=(v(i+1)-v(i))/delt;
end
derivative(L)=(-v(L))/delt;
end
