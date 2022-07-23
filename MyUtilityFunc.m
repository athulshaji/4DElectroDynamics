function MyUtilityFunc
clc;close all;clear variables;

IncFieldX1=load('IncFieldX1.txt');
IncFieldY1=load('IncFieldY1.txt');
IncFieldZ1=load('IncFieldZ1.txt');

IncFieldX2=load('IncFieldX2.txt');
IncFieldY2=load('IncFieldY2.txt');
IncFieldZ2=load('IncFieldZ2.txt');

IncFieldX3=load('IncFieldX3.txt');
IncFieldY3=load('IncFieldY3.txt');
IncFieldZ3=load('IncFieldZ3.txt');

IncFieldX4=load('IncFieldX4.txt');
IncFieldY4=load('IncFieldY4.txt');
IncFieldZ4=load('IncFieldZ4.txt');

FieldX1=load('FieldX1.txt');
FieldY1=load('FieldY1.txt');
FieldZ1=load('FieldZ1.txt');

FieldX2=load('FieldX2.txt');
FieldY2=load('FieldY2.txt');
FieldZ2=load('FieldZ2.txt');

FieldX3=load('FieldX3.txt');
FieldY3=load('FieldY3.txt');
FieldZ3=load('FieldZ3.txt');

FieldX4=load('FieldX4.txt');
FieldY4=load('FieldY4.txt');
FieldZ4=load('FieldZ4.txt');

delt=load('delt.txt');
T=load('T.txt');

Fs=1/delt;

Y = fft(FieldZ1,T);
FFTFieldZ1=Y;
L=T;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure,plot(f/1e9,((P1)));
title('FFT Ex');xlabel('GHz');ylabel('FFT mag');

Y = fft(FieldZ2,T);
FFTFieldZ2=Y;
L=T;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure,plot(f/1e9,((P1)));
title('FFT Ex');xlabel('GHz');ylabel('FFT mag');

Y=(FFTFieldZ2./FFTFieldZ1).^2;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
S21=P1;
f = Fs*(0:(L/2))/L;
figure,plot(f/1e9,((S21)));
title('S_{21}');xlabel('GHz');ylabel('FFT mag');
end