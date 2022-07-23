function HarmonicBalanceSolver_3_1
% Half wave rectifier with capacitor in shunt and resistor (shunt) as the
% load

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
DO_PLOT=0;
Rs=0.568; % Ohms
Isat=2.52e-9; % Amperes for Si diode
TempC=25; % degree Celcius
TempK=TempC+273.15; % Kelvin
Vth=1.38e-23*TempK/1.6e-19; %0.026 V for T=25degree Celcius;
eta=1.752;
K=25; % number of harmonics of interest k=0:K
SeriesResistorValue=50; % Ohms
ShuntResValue=500; % Ohms
ShuntCapValue=100e-9;
counter=1;
for dbm=-30:1:60
    clc
    dbm
    inputpowerinwatt=10^(dbm/10)/1000;
    SourceResValue=75;
    % Irms=sqrt(inputpowerinwatt/SourceResValue);
    % Vrms=Irms*SourceResValue;
    Vrms=sqrt(inputpowerinwatt*SourceResValue);
    Vpeak=Vrms*sqrt(2);
    T=1; % seconds
    omega0=50; % fundamental freq
    Fs=100*omega0;
    t=linspace(0,T,Fs*T);
    vsrc=Vpeak*sin(2*pi*omega0*t.'); % zeros(T,1);
    VSRC=fftshift(fft(vsrc));
    Y=(VSRC);
    L=length(vsrc);
    VSRC=Y/L;
    f=Fs*(-L/2:L/2-1)/L;
    vnonlinear=zeros(length(vsrc),1);% initial guess
    VNONLINEAR=fftshift(fft(vnonlinear));
    Y=VNONLINEAR;
    L=length(vnonlinear);
    lengthfreqvec=length(VNONLINEAR);
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
    
    VNONLINEARVEC=Y(freq_index)/L;
    VNONLINEARVEC=[VNONLINEARVEC;VNONLINEARVEC];
    
    x0=VNONLINEARVEC; % initial guess
    options = optimoptions('fsolve',...
        'Display','iter-detailed',...
        'Algorithm','levenberg-marquardt',...
        'FunctionTolerance',1e-9,...
        'MaxIterations',10000,...
        'MaxFunctionEvaluations',1000000,...
        'FiniteDifferenceType','central',...
        'PlotFcn',[]);
    [VOPTIM,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
    
    [VOPTIM(1:2*K+1) VOPTIM(2*K+2:end)]
    VNONLINEAR_NODE1=zeros(lengthfreqvec,1);
    VNONLINEAR_NODE2=zeros(lengthfreqvec,1);
    VNONLINEAR_NODE1(freq_index)=VOPTIM(1:2*K+1,1);
    VNONLINEAR_NODE2(freq_index)=VOPTIM(2*K+2:end,1);
    vnonlinear_node1=(ifft(ifftshift(VNONLINEAR_NODE1))*L);
    vnonlinear_node2=(ifft(ifftshift(VNONLINEAR_NODE2))*L);
    if(DO_PLOT)
        vdiode=vnonlinear_node1-vnonlinear_node2;
        for time=1:length(vdiode)
            if (vdiode(time)>=0)
                inonlinear_node1(time)=Isat*(exp(vdiode(time)/(eta*Vth))-1);
            elseif(vdiode(time)<0)
                inonlinear_node1(time)=-Isat;
            elseif(vdiode(time)==0)
                inonlinear_node1(time)=0;
            end
        end
        ilinear_node1=(vsrc-vnonlinear_node1)/SeriesResistorValue;
        figure,subplot 211
        plot(f,abs(VNONLINEAR_NODE1),'LineWidth',2);grid on; axis tight;ylabel('FFT amplitude');xlabel('freq');
        subplot 212
        plot(t,real(vnonlinear_node1),'LineWidth',2);grid on; axis tight;ylabel('Voltage amplitude');xlabel('Time(sec');
        figure,subplot 211
        plot(f,abs(VNONLINEAR_NODE2),'LineWidth',2);grid on; axis tight;ylabel('FFT amplitude');xlabel('freq');
        subplot 212
        plot(t,real(vnonlinear_node2),'LineWidth',2);grid on; axis tight;ylabel('Voltage amplitude');xlabel('Time(sec');
    end
    vmax2=max(real(vnonlinear_node2));
    vdc=vmax2/pi;
    idc=vdc/ShuntResValue;
    outputpower=idc^2*ShuntResValue;
    % rftodcconversionefficiency=((sum(abs(VOPTIM(2*K+2:end,1)).^2)/ShuntResValue)/sum(abs(VSRC).^2))*100
    Vdc=max(real(vnonlinear_node2))/pi;
    Idc=Vdc/ShuntResValue;
    percentagerftodcefficiency1(counter)=((Idc^2*ShuntResValue)/inputpowerinwatt)*100;
    % percentagerftodcefficiency2=((Idc^2*ShuntResValue)/inputpowerinwatt)*100
    clear negative_freq_index positive_freq_index freq_index
    counter=counter+1;
end
end


function F=CostFunction(VNONLINEARVEC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global lengthfreqvec;
global vsrc;
global freq_index;
global SeriesResistorValue ShuntCapValue ShuntResValue;
global Fs;
global Rs;
global Isat;
global Vth;
global eta;
global K;
global f;

L=length(vsrc);
VNONLINEAR_NODE1=zeros(lengthfreqvec,1);
VNONLINEAR_NODE2=zeros(lengthfreqvec,1);
VNONLINEAR_NODE1(freq_index)=VNONLINEARVEC(1:2*K+1,1);
VNONLINEAR_NODE2(freq_index)=VNONLINEARVEC(2*K+2:end,1);
vnonlinear_node1=(ifft(ifftshift(VNONLINEAR_NODE1))*L);
vnonlinear_node2=(ifft(ifftshift(VNONLINEAR_NODE2))*L);

% NODE1
ILINEAR_NODE1=fftshift(fft((vnonlinear_node1-vsrc)))/SeriesResistorValue; % current info from the linear part of the ckt
Y=ILINEAR_NODE1;
L=length(ILINEAR_NODE1);
ILINVEC1=(Y(freq_index)/L);
vdiode=vnonlinear_node1-vnonlinear_node2;
for time=1:length(vdiode)
    if (vdiode(time)>=0)
        inonlinear_node1(time)=Isat*(exp(vdiode(time)/(eta*Vth))-1);
    elseif(vdiode(time)<0)
        inonlinear_node1(time)=-Isat;
    elseif(vdiode(time)==0)
        inonlinear_node1(time)=0;
    end
end
INONLINEAR_NODE1=fftshift(fft(inonlinear_node1));
Y=INONLINEAR_NODE1;
L=length(INONLINEAR_NODE1);
INONLINVEC1=(Y(freq_index)/L).';

% NODE2
vdiode=vnonlinear_node1-vnonlinear_node2;
for time=1:length(vdiode)
    if (vdiode(time)>=0)
        inonlinear_node2(time)=-Isat*(exp(vdiode(time)/(eta*Vth))-1);
    elseif(vdiode(time)<0)
        inonlinear_node2(time)=Isat;
    elseif(vdiode(time)==0)
        inonlinear_node2(time)=0;
    end
end
INONLINEAR_NODE2=fftshift(fft(inonlinear_node2));
Y=INONLINEAR_NODE2;
L=length(INONLINEAR_NODE2);
INONLINVEC2=(Y(freq_index)/L).';
% Current from the other linear part of the circuit (RC section)
Z=@(f) ShuntResValue./(1+(1i*2*pi*f*ShuntResValue*ShuntCapValue));
IMPEDANCE=Z((f(freq_index))).';
ILINVEC2=(VNONLINEARVEC(2*K+2:end,1)./(IMPEDANCE));
F=[INONLINVEC1+ILINVEC1;...
    ILINVEC2+INONLINVEC2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
