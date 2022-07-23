function HarmonicBalanceSolver_3_2
% Half wave rectifier with capacitor in shunt and resistor (shunt) as the
% load.
% Working program!!!

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
dbm_array=[-25:5:20];
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
    Fs=1000*omega0;
    t=linspace(0,T,Fs*T);
    vsrc=Vpeak*cos(2*pi*omega0*t.'); % zeros(T,1);
    VSRC=fftshift(fft(vsrc));
    Y=(VSRC);
    L=length(vsrc);
    VSRC=Y/L;
    f=Fs*(-L/2:L/2-1)/L;
    vnonlinear=rand(length(vsrc),1);% initial guess
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
        'FunctionTolerance',1e-12,...
        'StepTolerance',1e-12,...
        'MaxIterations',100000,...
        'MaxFunctionEvaluations',1000000,...
        'FiniteDifferenceType','central',...
        'PlotFcn',[]);
    [VOPTIM,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
    exitflag
    [VOPTIM(1:2*K+1) VOPTIM(2*K+2:end)];
    VNONLINEAR_NODE1=zeros(lengthfreqvec,1);
    VNONLINEAR_NODE2=zeros(lengthfreqvec,1);
    VNONLINEAR_NODE1(freq_index)=VOPTIM(1:2*K+1,1);
    VNONLINEAR_NODE1_VEC=VOPTIM(1:2*K+1,1);
    VNONLINEAR_NODE2(freq_index)=VOPTIM(2*K+2:end,1);
    VNONLINEAR_NODE2_VEC=VOPTIM(2*K+2:end,1);
    vnonlinear_node1=(ifft(ifftshift(VNONLINEAR_NODE1))*L);
    vnonlinear_node2=(ifft(ifftshift(VNONLINEAR_NODE2))*L);
    if(REAL_DIODE_DATA)
        vdiode=real(vnonlinear_node1-vnonlinear_node2);
        for time=1:length(vdiode)
            if (vdiode(time)>0)
                inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                    real(Diode_Data.current),real(vdiode(time)),'spline');
                %         elseif(vdiode(time)<0 && vdiode(time)>=Bv)
                %             inonlinear_node1(time)=-Isat;
            elseif(vdiode(time)<0)
                inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                    real(Diode_Data.current),real(vdiode(time)),'spline');
            elseif(vdiode(time)==0)
                inonlinear_node1(time)=0;
            end
        end
        %     real_inonlinear_node1=interp1(real(Diode_Data.voltage),real(Diode_Data.current),real(vdiode),'spline');
        %     % imag_inonlinear_node1=interp1((Diode_Data.voltage),(Diode_Data.current),imag(vdiode),'spline');
        %     inonlinear_node1=real_inonlinear_node1;%+1i*imag_inonlinear_node1;
    end
    if(DIODE_EQN_BASED)
        vdiode=real(vnonlinear_node1-vnonlinear_node2);
        for time=1:length(vdiode)
            if (vdiode(time)>0)
                inonlinear_node1(time)=Isat*(exp(vdiode(time)/(eta*Vth))-1);
            elseif(vdiode(time)<0&&vdiode(time)>=Bv)
                inonlinear_node1(time)=-Isat;
            elseif(vdiode(time)<Bv)
                inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                    real(Diode_Data.current),real(vdiode(time)),'spline');
            elseif(vdiode(time)==0)
                inonlinear_node1(time)=0;
            end
        end
    end
%     ILINEAR_NODE1=(VSRC)/(SeriesResistorValue);
    ILINEAR_NODE2=VNONLINEAR_NODE2/ShuntResValue;
    ILINEAR_NODE1=fftshift(fft(inonlinear_node1))/length(inonlinear_node1);
    ILINEAR_NODE1_VEC=ILINEAR_NODE1(freq_index);
    ILINEAR_NODE2_VEC=ILINEAR_NODE2(freq_index);
    if(DO_PLOT)
        figure,subplot 211
        plot(f,abs(VNONLINEAR_NODE1),'LineWidth',2);grid on; axis tight;ylabel('FFT amplitude');xlabel('freq');
        subplot 212
        plot(t,real(vnonlinear_node1),'LineWidth',2);grid on; axis tight;ylabel('Voltage amplitude');xlabel('Time(sec');
        figure,subplot 211
        plot(f,abs(VNONLINEAR_NODE2),'LineWidth',2);grid on; axis tight;ylabel('FFT amplitude');xlabel('freq');
        subplot 212
        plot(t,real(vnonlinear_node2),'LineWidth',2);grid on; axis tight;ylabel('Voltage amplitude');xlabel('Time(sec');
    end
    percentage_RF2DC_conversion_efficiency(iter_number)=abs((real((VNONLINEAR_NODE2_VEC(K+1))*real(ILINEAR_NODE2_VEC(K+1))))/...
        abs(4*real((VNONLINEAR_NODE1_VEC(K+2))*conj(ILINEAR_NODE1_VEC(K+2))*0.5)))*100;
    disp(['dBm=',num2str(dbm),'; percentage efficiency=',num2str(percentage_RF2DC_conversion_efficiency),'%']);
    % percentagerftodcefficiency2=real(VNONLINEAR_NODE1()*conj(INONLINEAR_NODE1()*0.5)
    iter_number=iter_number+1;
    
end
toc;
% figure,plot(dbm_array,percentage_RF2DC_conversion_efficiency,'LineWidth',2);axis tight;grid on;
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
global Diode_Data
global REAL_DIODE_DATA DIODE_EQN_BASED
global Bv;

L=length(vsrc);
VNONLINEAR_NODE1=zeros(lengthfreqvec,1);
VNONLINEAR_NODE2=zeros(lengthfreqvec,1);
VNONLINEARVEC(K+1)=real(VNONLINEARVEC(K+1,1));
VNONLINEARVEC(2*K+1+K+1,1)=real(VNONLINEARVEC(2*K+1+K+1,1));
VNONLINEAR_NODE1(freq_index)=VNONLINEARVEC(1:2*K+1,1);
VNONLINEAR_NODE2(freq_index)=VNONLINEARVEC(2*K+2:end,1);
vnonlinear_node1=(ifft(ifftshift(VNONLINEAR_NODE1))*L);
vnonlinear_node2=(ifft(ifftshift(VNONLINEAR_NODE2))*L);

% NODE1
ILINEAR_NODE1=fftshift(fft((vnonlinear_node1-vsrc)))/SeriesResistorValue; % current info from the linear part of the ckt
Y=ILINEAR_NODE1;
L=length(ILINEAR_NODE1);
ILINVEC1=(Y(freq_index)/L);
if(REAL_DIODE_DATA)
    vdiode=real(vnonlinear_node1-vnonlinear_node2);
    for time=1:length(vdiode)
        if (vdiode(time)>0)
            inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
            %         elseif(vdiode(time)<0 && vdiode(time)>=Bv)
            %             inonlinear_node1(time)=-Isat;
        elseif(vdiode(time)<0)
            inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
        elseif(vdiode(time)==0)
            inonlinear_node1(time)=0;
        end
    end    
end
if(DIODE_EQN_BASED)
    vdiode=real(vnonlinear_node1-vnonlinear_node2);
    for time=1:length(vdiode)
        if (vdiode(time)>0)
            inonlinear_node1(time)=Isat*(exp(vdiode(time)/(eta*Vth))-1);
        elseif(vdiode(time)<0&&vdiode(time)>=Bv)
            inonlinear_node1(time)=-Isat;
        elseif(vdiode(time)<Bv)
            inonlinear_node1(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
        elseif(vdiode(time)==0)
            inonlinear_node1(time)=0;
        end
    end
end
INONLINEAR_NODE1=fftshift(fft(inonlinear_node1));
Y=INONLINEAR_NODE1;
L=length(INONLINEAR_NODE1);
INONLINVEC1=(Y(freq_index)/length(inonlinear_node1)).';

% NODE2
if(REAL_DIODE_DATA)
    vdiode=real(vnonlinear_node2)-real(vnonlinear_node1);
    for time=1:length(vdiode)
        if (vdiode(time)>0)
            inonlinear_node2(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
            %         elseif(vdiode(time)<0 && vdiode(time)>=Bv)
            %             inonlinear_node2(time)=-Isat;
        elseif(vdiode(time)<0)
            inonlinear_node2(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
        elseif(vdiode(time)==0)
            inonlinear_node2(time)=0;
        end
    end    
end
if(DIODE_EQN_BASED)
    vdiode=real(vnonlinear_node1-vnonlinear_node2);
    for time=1:length(vdiode)
        if (vdiode(time)>=0)
            inonlinear_node2(time)=-Isat*(exp(vdiode(time)/(eta*Vth))-1);
        elseif(vdiode(time)<0&&vdiode(time)>=Bv)
            inonlinear_node2(time)=Isat;
        elseif(vdiode(time)<Bv)
            inonlinear_node2(time)=interp1(real(Diode_Data.voltage),...
                real(Diode_Data.current),real(vdiode(time)),'spline');
        elseif(vdiode(time)==0)
            inonlinear_node2(time)=0;
        end
    end
end
INONLINEAR_NODE2=fftshift(fft(inonlinear_node2));
Y=INONLINEAR_NODE2;
L=length(INONLINEAR_NODE2);
INONLINVEC2=(Y(freq_index)/length(inonlinear_node1)).';
% Current from the other linear part of the circuit (RC section)
Z=@(f) ShuntResValue./(1+(1i*2*pi*f*ShuntResValue*ShuntCapValue));
IMPEDANCE=Z((f(freq_index))).';
ILINVEC2=(VNONLINEARVEC(2*K+2:end,1)./(IMPEDANCE));
F=[INONLINVEC1+ILINVEC1;...
    ILINVEC2+INONLINVEC2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

