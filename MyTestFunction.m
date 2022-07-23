function MyTestFunction
clc;close all;clear variables;
T=1; % seconds
omega0=50; % fundamental freq
Fs=100*omega0;
t=linspace(0,T,Fs);
vsrc_time=1*sin(2*pi*omega0*t.');
Y=fftshift(fft(vsrc_time));
L=length(vsrc_time);
P2=abs(Y/L);
f=Fs*(-L/2:L/2-1)/L;
signal=ifft(ifftshift(P2))*L;
figure
subplot 211
plot(f,P2);
subplot 212
plot(signal);
end