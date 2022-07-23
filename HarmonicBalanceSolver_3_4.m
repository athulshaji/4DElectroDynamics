function HarmonicBalanceSolver_3_4
% Developed by Athul Shaji, Microwave Lab, ECE, IISc Bangalore
% Half wave rectifier with capacitor in shunt and resistor (shunt) as the
% load.
% Developing a time-domain harmonic balance
% i1=linear current at node 1
% i1_=nonlinear current at node 1
% i2=linear current at node 2
% i2_=nonlinear current at node 2
% Trying to build a dynamic code which could be combined with FDTD

clc;clear all;close all;
format long
global lengthfreqvec;
global vsrc;
global positive_freq_index negative_freq_index freq_index;
global SeriesResistorValue ShuntCapValue ShuntResValue;
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

counter=1;
fid=fopen('SMS7630_Schottky_Diode.txt');
fgetl(fid);
while(~feof(fid))
    m=fgetl(fid);
    m=regexp(m,'	','split');
    m{1,1}=str2num(m{1,1});
    m{1,2}=str2num(m{1,2});
    Diode_Data.voltage(counter)=m{1};
    Diode_Data.current(counter)=m{2};
    counter=counter+1;
end
fclose(fid);
iter_number=1;
dbm_array=[0];
tic
for dbm=dbm_array
    tic
    clc
    dbm
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
    K=20; % number of harmonics of interest k=0:K
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
    Fs=1000*omega0;    
    t=linspace(0,T,Fs*T);
    vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);
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
    global time_index time K
    time=1:length(t);
    x0=[a1;a2]; % initial guess
    for time_index=1:length(t)
        options = optimoptions('fsolve',...
            'Display','iter-detailed',...
            'Algorithm','levenberg-marquardt',...
            'FunctionTolerance',1e-3,...
            'StepTolerance',1e-3,...
            'MaxIterations',100000,...
            'MaxFunctionEvaluations',1000000,...
            'FiniteDifferenceType','central',...
            'PlotFcn',[]);
        [a_optim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
        exitflag
    end
end
a1=a_optim(1:2*K+1);
a2=a_optim(2*K+2:end);
for i=1:2*K+1
    v1(i,:)=a1(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));    
end
for i=1:2*K+1    
    v2(i,:)=a2(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t));
end
v1=sum(v1);v2=sum(v2);
figure('Color','w');
subplot 211
plot(t,v1,'ok');axis tight;grid on;
subplot 212
plot(t,v2,'-ok');axis tight;grid on;
end

function F=CostFunction(a)
global time K
global lengthfreqvec;
global vsrc;
global positive_freq_index negative_freq_index freq_index;
global SeriesResistorValue ShuntCapValue ShuntResValue;
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
global time_index

a1=a(1:2*K+1);
a2=a(2*K+2:end);
% Node 1
for i=1:2*K+1
    v1(i,:)=a1(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t(1:time_index)));    
end
for i=1:2*K+1    
    v2(i,:)=a2(i).*exp(1i*2*pi.*(f(freq_index(i))).*(t(1:time_index)));
end
v1=real(sum(v1));
v2=real(sum(v2));
i1=(v1-vsrc(1:time_index).')/SeriesResistorValue;
i1_temp=zeros(1,length(t));
for i=1:length(v2)
   i1_temp(i)=i1(i);    
end
i1=i1_temp;
vdiode=v1-v2;
for T=time_index
    if (vdiode(T)>0)
        i1_(T)=Isat*(exp(vdiode(T)/(eta*Vth))-1);
    elseif(vdiode(T)<0&&vdiode(T)>=Bv)
        i1_(T)=-Isat;
    elseif(vdiode(T)<Bv)
        i1_(T)=interp1(real(Diode_Data.voltage),...
            real(Diode_Data.current),real(vdiode(T)),'spline');
    elseif(vdiode(T)==0)
        i1_(T)=0;
    end
end
i1_temp=zeros(1,length(t));
for i=1:length(i1_)
   i1_temp(i)=i1_(i);    
end
i1_=i1_temp;
% Node 2
vdiode=v1-v2;
for T=time_index
    if (vdiode(T)>=0)
        i2_(T)=Isat*(exp(vdiode(T)/(eta*Vth))-1);
    elseif(vdiode(T)<0&&vdiode(T)>=Bv)
        i2_(T)=-Isat;
    elseif(vdiode(T)<Bv)
        i2_(T)=interp1(real(Diode_Data.voltage),...
            real(Diode_Data.current),real(vdiode(T)),'spline');
    elseif(vdiode(T)==0)
        i2_(T)=0;
    end
end
i2_temp=zeros(1,length(t));
for i=1:length(i2_)
   i2_temp(i)=i2_(i);    
end
i2_=i2_temp;
v2_temp=zeros(1,length(t));
for i=1:length(v2)
   v2_temp(i)=v2(i);    
end
v2=v2_temp;
V2=fftshift(fft(v2))/length(v2);
V2=V2(freq_index);
Z=@(f) ShuntResValue./(1+(1i*2*pi*f*ShuntResValue*ShuntCapValue));
IMPEDANCE=Z((f(freq_index)));
I2=(V2)./IMPEDANCE;
I22=zeros(1,length(vsrc));
I22(freq_index)=I2;
i2=ifft(ifftshift((I22)))*length(vsrc);
F=[i2+i2_;...
    i1+i1_];
clear v1 v2
end


