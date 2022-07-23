clc;close all;clear all;
format long
Fs=50000;
t=0:(1/Fs):1;
x=sin(2*pi*1000*t);
L=length(x);
f=Fs*(-L/2:(L/2-1))/L;
X=fftshift(fft(x))/L;
figure,
subplot 211 
plot(t,x),axis tight;grid on;
subplot 212
stem(f,abs(X));grid on;axis tight;
figure,stft(x,Fs);

