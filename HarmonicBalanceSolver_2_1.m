function HarmonicBalanceSolver_2_1
% voltage source in shunt and the Resistor in series followed by a shunt
% diode (1N4001--> Si p-n junction diode). An inductor in series with the capacitor is placed in shunt
% with the diode. This is the modified code for doing statistical studies

clc;clear all;close all;
format long
for numHarmonics=1:50
    clc
    numHarmonics
    global lengthfreqvec;
    global vsrc_time;
    global positive_freq_index negative_freq_index freq_index;
    global R Cval Lval;
    global Fs;
    global Rs;
    global Isat;
    global Vth;
    global eta;
    global K;
    global f;
    DO_PLOT=0;
    Rs=0.568; % Ohms
    Isat=2.52e-9; % Amperes for Si diode
    TempC=25; % degree Celcius
    TempK=TempC+273.15; % Kelvin
    Vth=1.38e-23*TempK/1.6e-19; %0.026 V for T=25degree Celcius;
    eta=1.752;
    K=numHarmonics; % number of harmonics of interest k=0:K
    R=50; % Ohms
    Lval=4.7e-9;
    Cval=100e-9;
    T=1; % seconds
    omega0=50; % fundamental freq
    Fs=1000*omega0;
    t=linspace(0,T,Fs*T);
    vsrc_time=1*sin(2*pi*omega0*t.'); % zeros(T,1);
    Y=fftshift(fft(vsrc_time));
    L=length(vsrc_time);
    P2=abs(Y/L);
    f=Fs*(-L/2:L/2-1)/L;
    % P1=P2(1:L/2+1);
    % P1(2:end-1)=2*P1(2:end-1);
    % f=Fs*(0:(L/2))/L;
    vnonlinear_time=zeros(length(vsrc_time),1);% initial guess
    vnonlinear_freq=fftshift(fft(vnonlinear_time));
    lengthfreqvec=length(vnonlinear_freq);
    Y=vnonlinear_freq;
    L=length(vnonlinear_time);
    P2=abs(Y/L);
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
    vnonlinearvec=Y(freq_index)/L;
    x0=vnonlinearvec; % initial guess
    options = optimoptions('fsolve',...
        'Display','iter-detailed',...
        'Algorithm','levenberg-marquardt',...
        'FunctionTolerance',1e-9,...
        'MaxIterations',1000,...
        'MaxFunctionEvaluations',100000,...
        'FiniteDifferenceType','central');
    % voptim=fminunc(@CostFunction,x0,options);
    [voptim{numHarmonics},fval,exitflag,output]=fsolve(@CostFunction,x0,options);
    %     filename=['voptim',num2str(numHarmonics),'.txt'];
    
    if(DO_PLOT)
        % voptim=fftshift(voptim);
        Z=@(f) (1i*2*pi*f).*(Lval-(1./(Cval.*((2*pi*f).^2))));
        IMPEDANCE=Z(f(freq_index)).';
        IMPEDANCE(isnan(IMPEDANCE))=-1i*1e20;
        ilinvec2=((voptim)./(IMPEDANCE));
        ilinear_freq_2=zeros(lengthfreqvec,1);
        ilinear_freq_2(freq_index)=ilinvec2;
        voptim
        voptim_freq=zeros(lengthfreqvec,1);
        voptim_freq(freq_index)=voptim;
        [pks,locs] = findpeaks(abs(voptim_freq));
        figure
        subplot 211
        plot(f,abs(voptim_freq),f(locs),pks,'-or','LineWidth',2);title(['K=',num2str(K)]);grid on
        xlim([-K*omega0 K*omega0]);
        xlabel('Freq(Hz)');ylabel('FFT Mag');
        subplot 212
        plot(t,real(ifft(ifftshift(voptim_freq))*L),'LineWidth',2);grid on;
        figure,
        subplot 211
        plot(f,abs(ilinear_freq_2),'LineWidth',2);
        subplot 212
        plot(t,real(ifft(ifftshift(ilinear_freq_2))*L),'LineWidth',2),grid on;
        
        vnonlinear_freq=zeros(lengthfreqvec,1);
        % vnonlinearvec(K+1)=real(vnonlinearvec(K+1))+1i*0;
        vnonlinear_freq(freq_index)=voptim;
        vnonlinear_time=(ifft(ifftshift(vnonlinear_freq))*L);
        figure,plot(t,real(vnonlinear_time-vsrc_time)/R,'LineWidth',2);grid on;hold on
        
        
        vnonlinear_time=(ifft(ifftshift(voptim_freq))*L);
        % x0=0;
        % global vnonlinear_time time
        % for time=1:length(t)
        %     [idiode_time,fval2,exitflag2,output2]=fsolve(@DiodeCurrent,x0,options);
        %     diode_current(time)=idiode_time;
        % end
        % figure,plot(t,diode_current,'LineWidth',2),grid on;
        % legend('current through R','current through diode');
        for time=1:length(vnonlinear_time)
            if (vnonlinear_time(time)>=0)
                inonlinear_time(time)=Isat*(exp(vnonlinear_time(time)/(eta*Vth))-1);
            elseif(vnonlinear_time(time)<0)
                inonlinear_time(time)=-Isat;
            elseif(vnonlinear_time(time)==0)
                inonlinear_time(time)=0;
            end
        end
        plot(t,real(inonlinear_time),'LineWidth',2);
        legend('current through R','current through diode');
    end
    exitflag
end
for i=5:numHarmonics
    absvoptim{i}=abs(voptim{i});
    dataDC(i)=absvoptim{i}(i+1);% mag of DC value
    dataFUN(i)=absvoptim{i}(i+2);% mag of Fundamental freq
    data2Har(i)=absvoptim{i}(i+3);% mag of 2nd harmonic freq
    data3Har(i)=absvoptim{i}(i+4);% mag of 3rd harmonic freq
    data4Har(i)=absvoptim{i}(i+5);% mag of 4th harmonic freq
    data5Har(i)=absvoptim{i}(i+6);% mag of 5th harmonic freq
end
end

% Cost function for the HB problem
function F=CostFunction(vnonlinearvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global lengthfreqvec;
global vsrc_time;
global positive_freq_index negative_freq_index freq_index;
global R Cval Lval;
global Fs;
global Isat;
global Vth;
global eta;
global K;
global f;
global time vnonlinear_time;
L=length(vsrc_time);
vnonlinear_freq=zeros(lengthfreqvec,1);
vnonlinear_freq(freq_index)=vnonlinearvec;
vnonlinear_time=(ifft(ifftshift(vnonlinear_freq))*L);
ilinear_freq=fftshift(fft((vnonlinear_time-vsrc_time))/R); % current info from the linear part of the ckt
Y=ilinear_freq;
L=length(vsrc_time);
P2=abs(Y/L);
ilinvec=(Y(freq_index)/L);
% Current from the other linear part of the circuit (LC section)
Z=@(f) (1i*2*pi*f).*(Lval-(1./(Cval.*((2*pi*f).^2))));
IMPEDANCE=Z((f(freq_index)))';
IMPEDANCE(isnan(IMPEDANCE))=-1i*1e20; % replacing infinity with a large valuefor imaginary part
ilinvec2=(vnonlinearvec./(IMPEDANCE));
options = optimoptions('fsolve',...
    'Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-2,...
    'MaxIterations',1000,...
    'MaxFunctionEvaluations',100000);
% for time=1:L
%     x0=0;
%     [idiode_time,fval2,exitflag2,output2]=fsolve(@DiodeCurrent,x0,options);
%     diode_current(time)=idiode_time;
% end
% inonlinear_freq=fftshift(fft(diode_current));
for time=1:length(vnonlinear_time)
    if (vnonlinear_time(time)>=0)
        inonlinear_time(time)=Isat*(exp(vnonlinear_time(time)/(eta*Vth))-1);
    elseif(vnonlinear_time(time)<0)
        inonlinear_time(time)=-Isat;
    elseif(vnonlinear_time(time)==0)
        inonlinear_time(time)=0;
    end
end
inonlinear_freq=fftshift(fft(inonlinear_time));
Y=inonlinear_freq;
L=length(vsrc_time);
inonlinvec=(Y(freq_index)/L).';
F=(inonlinvec+ilinvec+ilinvec2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Cost function for the diode current
function error=DiodeCurrent(diode_current)
global lengthfreqvec;
global vsrc_time;
global positive_freq_index negative_freq_index freq_index;
global R Cval Lval;
global Fs;
global Rs
global Isat;
global Vth;
global eta;
global K;
global f;
global vnonlinear_time time;
error=Isat*(exp((real(vnonlinear_time(time))-diode_current*Rs)/(eta*Vth))-1)-diode_current;
end

