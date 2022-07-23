function HarmonicBalanceSolver_1
% voltage source in shunt and the Resistor in series followed by a shunt
% diode (1N4001--> Si p-n junction diode)

clc;clear all;close all;
format long
global lengthfreqvec;
global vsrc_time;
global positive_freq_index negative_freq_index freq_index;
global R;
global Fs;
global Isat;
global Vth;
global eta;
global K;
Isat=1e-14; % Amperes for Si diode
TempC=25; % degree Celcius
TempK=TempC+273.15; % Kelvin
Vth=1.38e-23*TempK/1.6e-19; %0.026 V for T=25degree Celcius;
eta=1.05;
K=20; % number of frequencies of interest
R=50; % Ohms
T=1; % seconds
omega0=50; % fundamental freq
Fs=100*omega0;
t=linspace(0,T,Fs*T);
vsrc_time=2*sin(2*pi*omega0*t.'); % zeros(T,1);
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
    'FunctionTolerance',1e-3,...
    'MaxIterations',1000);
% voptim=fminunc(@CostFunction,x0,options);
[voptim,fval,exitflag,output]=fsolve(@CostFunction,x0,options);
% voptim(K+1)=real(voptim(K+1))+1i*0;
voptim%=sort(voptim,'descend')
voptim_freq=zeros(lengthfreqvec,1);
voptim_freq(freq_index)=voptim;
[pks,locs] = findpeaks(abs(voptim_freq));
figure
subplot 211
plot(f,abs(voptim_freq),f(locs),pks,'-or','LineWidth',2);title(['K=',num2str(K)]);grid on
xlim([0 K*omega0]);
xlabel('Freq(Hz)');ylabel('FFT Mag');
subplot 212
plot(t,real(ifft(ifftshift(voptim_freq))*L),'LineWidth',2);grid on;
v1=real(ifft(ifftshift(voptim_freq))*L);
figure,
subplot 211
plot(t,vsrc_time,t,v1,'LineWidth',2);grid on
subplot 212
plot(t,(vsrc_time-v1)/R,'LineWidth',2);grid on;
exitflag
end

% Cost function for the HB problem
function error=CostFunction(vnonlinearvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global lengthfreqvec;
global vsrc_time;
global positive_freq_index negative_freq_index freq_index;
global R;
global Fs;
global Isat;
global Vth;
global eta;
global K;
L=length(vsrc_time);
vnonlinearvec=complex(vnonlinearvec);
vnonlinear_freq=zeros(lengthfreqvec,1);
% vnonlinearvec(K+1)=real(vnonlinearvec(K+1))+1i*0;
vnonlinear_freq(freq_index)=vnonlinearvec;
vnonlinear_time=(ifft(ifftshift(vnonlinear_freq))*L);
ilinear_freq=fftshift(fft((vnonlinear_time-vsrc_time)/R)); % current info from the linear part of the ckt
Y=ilinear_freq;
L=length(vsrc_time);
P2=abs(Y/L);
% P1=P2(1:L/2+1);
% P1(2:end-1)=2*P1(2:end-1);
% f=Fs*(0:(L/2))/L;
ilinvec=(Y(freq_index)/L);
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
L=length(inonlinear_time);
P2=abs(Y/L);
% P1=P2(1:L/2+1);
% P1(2:end-1)=2*P1(2:end-1);
% f=Fs*(0:(L/2))/L;
inonlinvec=(Y(freq_index)/L);
error=((inonlinvec+ilinvec)); % we have to get this as zero. For KCL total current at a node is ZERO
% we need to minimize errorfunc to get the correct value of currents and
% voltages at the nodes
% error=sum(abs((error)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

